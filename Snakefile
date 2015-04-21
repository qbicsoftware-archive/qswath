import os
import subprocess
from os.path import exists as pexists
from os.path import join as pjoin
import hashlib


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
INI_PATH = os.path.join(SNAKEDIR, 'inis')

DATA = config['data']
RESULT = config['result']
LOGS = config['logs']
REF = config['ref']

try:
    path = subprocess.check_output(["which", "IDMerger"]).decode()
    OPENMS_BIN = os.path.dirname(path)
except subprocess.CalledProcessError:
    OPENMS_BIN = "/usr/bin"

class OpenMS:
    def __init__(self, bin_path, ini_dir, log_dir):
        self._bin_path = bin_path
        self._tools = os.listdir(bin_path)
        self._ini_path = ini_dir
        if not pexists(log_dir):
            os.mkdir(log_dir)
        self._log_dir = log_dir

    def __getattr__(self, name):
        if name not in self._tools:
            raise ValueError("OpenMS tool not found: %s" % name)

        def wrapper(input, output, logfile=None, extra_args=None, ini=None):
            if output is not None and not isinstance(output, list):
                if not pexists(os.path.dirname(str(output))):
                    os.mkdir(os.path.dirname(str(output)))
            elif output is not None:
                for out in output:
                    if not pexists(os.path.dirname(out)):
                        os.mkdir(os.path.dirname(out))
            if extra_args is None:
                extra_args = []
            if logfile is None:
                if output is None:
                    identifier = str(uuid.uuid4())
                else:
                    identifier = os.path.basename(str(output))
                logfile = '{}_{}'.format(name, identifier)

            command = [pjoin(self._bin_path, name)]
            if input is not None:
                if isinstance(input, list):
                    command += ['-in'] + list(input)
                else:
                    command += ['-in'] + [str(input)]
            if output is not None:
                if isinstance(output, list):
                    command += ['-out'] + list(output)
                else:
                    command += ['-out'] + [str(output)]
            if ini is not None:
                if isinstance(ini, list):
                    if not len(ini) == 1:
                        raise ValueError("Invalid params.")
                    ini = ini[0]
                command += ['-ini', str(ini).split(',', 1)[1]]
            command += extra_args

            log_std = pjoin(self._log_dir, logfile + '.out')
            log_err = pjoin(self._log_dir, logfile + '.err')
            with open(log_std, 'w') as out:
                with open(log_err, 'w') as err:
                    subprocess.check_call(command, stdout=out, stderr=err)

        return wrapper

openms = OpenMS(OPENMS_BIN, INI_PATH, LOGS)

# store the content of the ini file, so that snakemake will run
# rules agrain if the parameters inside the file change.
def params(name):
    path = pjoin(INI_PATH, name + '.ini')
    try:
        with open(path, 'rb') as f:
            # TODO replace makes sure that there are no wildcards
            # in the file content. This is not exactly clean ;-)
            return "{},{}".format(hashlib.sha256(f.read()).hexdigest(), path)
    except FileNotFoundError as e:
        raise ValueError("ini file '%s' not found" % path) from e


def make_ini_diff():
    orig_ini = pjoin(R_HOME, '..', 'inis')
    ini_diff = subprocess.Popen(
        ['diff', '-u', '-w', orig_ini, INI_PATH],
        stdout=subprocess.PIPE
    )
    ini_diff.wait()
    return ini_diff.communicate()[0].decode()


DDA_INPUT = os.path.join(config['data'],
                         config['params']['mzml_dda'] + ".mzML")
DIA_INPUT = config['params']['mzml_dia']
WINDOWS = os.path.join(config['data'],
                       config['params']['windows'])


for name in DIA_INPUT:
    file = os.path.join(config['data'], name + ".mzML")
    if not os.path.exists(file):
        raise ValueError("Could not find file %s" % file)


if not os.path.exists(DDA_INPUT):
    file = os.path.join(config['data'], name + ".mzML")
    raise ValueError("Could not find file %s" % DDA_INPUT)


if not os.path.exists(WINDOWS):
    raise ValueError("Could not find file %s" % WINDOWS)


rule all:
    input: expand("{result}/peptides.csv", name=DDA_INPUT, result=RESULT),
           expand("{result}/report.html", result=RESULT)


rule report:
    input: expand("{result}/peptides.csv", result=RESULT)
    output: expand("{result}/report.html", result=RESULT)
    run:
        with open(str(output), 'w') as f:
            json.dump(f, {"config": config, "ini_diff": make_ini_diff(),
                          "version": VERSION})


rule decoy:
    input: [pjoin(config['ref'], fasta) for fasta in config['params']['fasta']]
    output: "Decoy/database.fasta"
    shell:
        "cat {input} | decoyFastaGenerator.pl - DECOY_ - > {output}"


rule xinteract:
    input: db="Decoy/database.fasta", mzml=DDA_INPUT
    output: "Search/peptides.pep.xml"
    shell: "xinteract -N{output} -i -dDECOY_ -OAdtP {input}"


rule spectrast:
    input: "Search/peptides.pep.xml"
    output: "Filter/peptides.pep.xml"
    shell: "spectrast -c_BIN! -cf'Protein!~DECOY' -cN{output} -cP0.001 {input}"


rule spectrast_irt:
    input: "Filter/peptides.pep.xml"
    output: "IRT/peptides.splib"
    shell: "spectrast2spectrast_irt.py -i {input} -o {output}"


rule consensus:
    input: "IRT/peptides.splib"
    output: "Consensus/peptides.splib"
    shell: "spectrast -c_BIN! -cAC -cN{output} {input}"


rule to_tsv:
    input: "Consensus/peptides.splib"
    output: "TSV/library.csv"
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


rule ExtractChromatogram:
    input: mzml=os.path.join(config['data'], "{name}.mzML"), \
           lib="library.TraML"
    output: "ExtractChromatogram/{name}.mzML"
    params: params("OpenSwathChromatogramExtractor")
    run:
        extra_args = ['-tr', input['lib'], '-is_swath']
        openms.OpenSwathChromatogramExtractor(
            input['mzml'], output, ini=params, extra_args=extra_args
        )


rule RTNormalize:
    input: expand("ExtractChromatogram/{name}.mzML", name=DIA_INPUT), \
           lib="library.TraML"
    output: "RTNormalized/trafo.trafoXML"
    params: params('OpenSwathRTNormalizer')
    run:
        merged = "RTNormalized/_merged.mzML"
        openms.FileMerger(input, merged)

        extra_args = ['-tr', input['lib']]
        openms.OpenSwathRTNormalizer(
            merged, output, ini=params, extra_args=extra_args
        )
        os.remove(merged)


rule ExtractChromatogramNorm:
    input: mzml=os.path.join(config['data'], "{name}.mzML"), \
           lib="library.TraML", \
           trafo="RTNormalized/trafo.trafoXML"
    output: "ExtractChromatogramNorm/{name}.mzML"
    params: params("OpenSwathChromatogramExtractor")
    run:
        extra_args = ['-tr', input['lib'], '-is_swath']
        extra_args += ['-rt_norm', input['trafo']]
        openms.OpenSwathChromatogramExtractor(
            input['mzml'], output, ini=params, extra_args=extra_args
        )


rule OpenSwathAnalyzer:
    input: normed="ExtractChromatogramNorm/{name}.mzML", \
           mzml=os.path.join(config['data'], "{name}.mzML")
    output: "Analysed/{name}.feature.XML"
    params: params("OpenSwathAnalyzer")
    run:
        extra_args = ['-swath_files', input['mzml']]
        openms.OpenSwathAnalyzer(
            input['normed'], output, ini=params, extra_args=extra_args
        )


rule ToTSV:
    input: features=[os.path.join("Analysed", name + ".feature.XML") \
                     for name in DIA_INPUT], \
           mzml=[os.path.join("ExtractChromatogramNorm", name + '.mzML') \
                 for name in DIA_INPUT]
    output: os.path.join(config['result'], "peptides.csv")
    run:
        merged_feature = "ToTSV/merged.featureXML"
        openms.FileMerger(input['features'], merged_feature)
        openms.OpenSwathFeatureXMLtoTSV(
            merged_feature, output
        )
