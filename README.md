# VLBI-utilities
Utilities that can be used at stations to help with VLBI observations. For example:


Generation of antab files for continuous calibration (80 Hz). Usage:

```bash
antabfs.py fs_log_file
```


Checking if schedules at the station are updated. The script logs in vlbeer and checks the VEX files at the station against the version in the FTP site and shows the results. If the file is outdated it warns the user. Usage:

```bash
checkvexversions.py  monthyear stationCode     (Example: checkvexversions.py oct15 Ys)
```

Checking the quality of amplitude calibration after correlation and analysis. The data are downloaded from JIVE web page: http://archive.jive.nl/exp/. Usage:

```bash
getampcal.py year month stationCode    (Example: get_ampcal.py 2015 10 YS)
```

Script for transferring GPS files to Vlbeer?
