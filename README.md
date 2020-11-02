# VLBI-utilities
Utilities that can be used at stations to help with VLBI observations. For example:


## Generation of ANTAB files

Generation of antab files for continuous calibration (80 Hz). Usage:

```bash
antabfs.py [-f rxg_files_list] fs_log_file
```
All RXG files are supposed to be under /usr2/control/rxg_files/. If the -f option is not given, the script
will search there for a valid RXG file. Valid files are those which define a frequency range that contains 
the observed setup in the log file AND match the station code in the log file name. To do so it must be
named with the station code as Sc, e.g.:
calYsQ.rxg

Once the ANTAB files are generated at the station, please upload them to the `vlbeer.ira.inaf.it` ftp server. Place the files inside the `vlbi_arch/mmmYY/` folder, where `mmmYY` is the *date of the observation* (e.g. `mar20` if the experiment was observed in March 2020).


## VEX files up-to-date?

Checking if schedules at the station are updated. The script logs in vlbeer and checks the VEX files at the station against the version in the FTP site and shows the results. If the file is outdated it warns the user. Usage:

```bash
checkvexversions.py  monthyear stationCode     (Example: checkvexversions.py oct15 Ys)
```


## Amplitude Calibration Quality

Checking the quality of amplitude calibration after correlation and analysis. The data are downloaded from JIVE web page: http://archive.jive.nl/exp/. Usage:

```bash
getampcal.py year month stationCode    (Example: get_ampcal.py 2015 10 YS)
```

Script for transferring GPS files to Vlbeer?
