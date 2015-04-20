A workflow that creates spectral libraries for OpenSwath and uses it
to analyse data.

This is very much *not* finished!

Dependencies:

- OpenMS >= 2.0
- msproteomicstools
- TPP
- snakemake

All tools must be in the PATH.

The config file (`config.json`) should look like this:

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
	"fasta": ["database.fasta"]
    }
}
```
