# qswath

A workflow that creates spectral libraries for OpenSwath and uses it
to analyse data.

This is very much *not* finished!

## Quickstart

Install anaconda.

Run the following commands:

```
pip install git+https://github.com/qbicsoftware/qproject
cd your/workdir
qproject create -t . -w github:qbicsoftware/qswath
cp src/inis/* etc
bash src/build
vim src/config.json  # add a section "params" with the parameters you need
cp your/input/data data
cp your/fasta ref
qproject run -t .
```

## Config file

The config file (`config.json`) should look something like this:

```json
{
    "base": "/home/adr/git/QBiC/qswath/testproj",
    "data": "/home/adr/git/QBiC/qswath/testproj/data",
    "ref": "/home/adr/git/QBiC/qswath/testproj/ref",
    "src": "/home/adr/git/QBiC/qswath/testproj/src/qswath",
    "var": "/home/adr/git/QBiC/qswath/testproj/var/qswath",
    "result": "/home/adr/git/QBiC/qswath/testproj/result/qswath",
    "run": "/home/adr/git/QBiC/qswath/testproj/run/qswath",
    "inis": "/home/adr/git/QBiC/qswath/testproj/inis/qswath",
    "logs": "/home/adr/git/QBiC/qswath/testproj/logs/qswath",
    "params": {
        "mzml_dda": "lib",
        "mzml_dia": ["data", "data2"],
        "windows": "windows.csv",
        "fasta": ["database.fasta", "crap.fasta"]
    }
}
```
