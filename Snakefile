import os
import subprocess
from os.path import exists as pexists
from os.path import join as pjoin
import hashlib
from glob import glob


configfile: "config.json"

workdir: config['var']
SNAKEDIR = config['src']

try:
    VERSION = subprocess.check_output(
        ['git', 'describe', '--tags', '--always', '--dirty'],
        cwd=SNAKEDIR
    ).decode().strip()
except subprocess.CalledProcessError:
    VERSION = 'unknown'

R_HOME = os.path.join(SNAKEDIR, 'r_scripts')
DEFAULT_INIS = os.path.join(SNAKEDIR, 'inis')

DATA = config['data']
RESULT = config['result']
LOGS = config['logs']
REF = config['ref']
ETC = config['etc']


def etc(path):
    return os.path.join(ETC, path)


def logs(path):
    return os.path.join(LOGS, path)


def data(path):
    return os.path.join(DATA, path)


def ref(path):
    return os.path.join(REF, path)


def result(path):
    return os.path.join(RESULT, path)

params_path = etc('global.json')
if os.path.exists(params_path):
    with open(params_path) as f:
        config['params'] = json.load(f)
else:
    config['params'] = {}

if 'dda_data' not in config['params']:
    names = [file[:-len('.mzML')] for file in glob(data('*.mzML'))]
    names = [os.path.split(name)[-1] for name in names]
    config['params']['dda_data'] = names
if 'fasta' not in config['params']:
    config['params']['fasta'] = glob(ref('*.fasta'))

DDA_INPUT = config['params']['dda_data']
WINDOWS = ref('windows.txt')

for name in DDA_INPUT:
    if not os.path.exists(data(name) + '.mzML'):
        raise ValueError("Could not find file %s" % file)


if not os.path.exists(WINDOWS):
    raise ValueError("Could not find file %s" % WINDOWS)


rule all:
    input: expand("{result}/library.csv", name=DDA_INPUT, result=RESULT),
           expand("{result}/report.html", result=RESULT)


rule report:
    output: expand("{result}/report.html", result=RESULT)
    run:
        with open(str(output), 'w') as f:
            json.dump({"config": config, "version": VERSION}, f)


rule decoy:
    input: config['params']['fasta']
    output: "Decoy/database.fasta"
    shell:
        "cat {input} | decoyFastaGenerator.pl - DECOY_ - > {output}"


rule ConvertMZ5:
    input: data('{name}.mzML')
    output: "MZXML/{name}.mzXML"
    run: "msconvert {input} --mzXML -o MZ5"


rule comet:
    input:
        mzxml="MZXML/{name}.mzXML",
        fasta="Decoy/database.fasta"
    output: "Search/comet_{name}.pep.xml"
    run:
        params = etc("comet.params")
        name = wildcards['name']
        with open(os.path.join(LOGS, "comet_" + name), "w") as f:
            subprocess.check_call(
                [
                    "comet",
                    "-P" + params,
                    "-D" + str(input['fasta']),
                    "-N{}".format(str(output)[:-len('.pep.xml')]),
                    str(input['mzxml'])
                ],
                stdout=f,
                stderr=f,
            )


rule xinteract:
    input: expand("Search/comet_{name}.pep.xml", name=DDA_INPUT)
    output: "PeptideProphet/peptides.pep.xml"
    shell: "xinteract -N{output} -i -dDECOY_ -OAdtP {input}"



rule spectrast:
    input: "PeptideProphet/peptides.pep.xml"
    output: "Spectrast/peptides.splib"
    shell:
        "spectrast -c_BIN! -cf'Protein!~DECOY' -cNtmp.pep.xml -cP0.001 {input}"
        "spectrast2spectrast_irt.py -i tmp.pep.xml -o {output}"


rule consensus:
    input: "Spectrast/peptides.splib"
    output: "Consensus/peptides.splib"
    shell: "spectrast -c_BIN! -cAC -cN{output} {input}"


rule to_tsv:
    input: "Consensus/peptides.splib"
    output: result("library.csv")
    shell:
        "spectrast2tsv.py -n 6 -o 6 -s y,b -k openswath " +
            "-w " + WINDOWS + " -a {output} {input}"


rule to_TraML:
    input: "TSV/library.csv"
    output: "TraML/library.TraML"
    shell:
        "ConvertTSVToTraML -in {input} -out {output}"


rule swath_decoys:
    input: "TraML/library.TraML"
    output: "library.TraML"
    shell:
        "OpenSwathDecoyGenerator -in {input} -out {output} -method shuffle "
        "-append -exclude_similar"
