"""External lookup file for sensors, header information.
This file may need to be edited at sea if a sensor is attached to a package that
this software has not encountered before. In that case, follow the steps below:

1. On the SBE acquisition machine, open the .XMLCON for a cast.
2. Look down the .XMLCON file until you find a line with:
        <Sensor index="xx" SensorID="yy" >
   Where "xx" and "yy" are numbers, and the following lines are the sensor to edit in.
3. In the short_lookup table below, copy and paste a line, then edit the line:
        'SensorID':{'short_name': 'mm', 'long_name':'ll', 'units': 'uu', 'type': 'tt'}
    Where mm is a science defined short name, ll is the name of the sensor,
    uu is the unit used by the sensor, and tt is the data, computer format of the sensor,
    which will typically be float64. (Any sensor that puts out a decimal number is float64)
4. If adding at the end of the short_lookup table, remember to add a comma ',' to the
    previous line, if adding at the top or middle add a comma to the inserted line.
5. Save the file, and attempt to rerun the program again.

Questions for the future:
-This is SBE focused right now. What do we do for potential namespace conflicts with
other CTDs/sensors from other companies (RBR, others)?
-CTDs send up signals and could have varying levels of precision. Do we put a default
truncation value in a lookup table, or make sure the corresponding equation has it built in
(which again with SBE CTDs are dependent on how the CTD logs it, unless it has a datalogger),
or a third option?

Joseph Gum
November 10, 2016
"""

#lookup table for sensor data
###DOUBLE CHECK TYPE IS CORRECT###
short_lookup = {
    '55':{'short_name': 't', 'long_name':'SBE 3+ Temperature', 'units': 'C', 'type': 'float64'},
    '45':{'short_name': 'p', 'long_name':'SBE 9+ Pressure', 'units': 'dbar', 'type': 'float64'},
    '3':{'short_name': 'c', 'long_name':'SBE 4 Conductivity', 'units': 'S/m', 'type':'float64'},
    '38':{'short_name': 'o', 'long_name':'SBE 43 Oxygen', 'units': 'ml/l', 'type':'float64'},
    '11':{'short_name': 'fluoro', 'long_name':'Seapoint Fluorometer', 'units': 'ug/l', 'type':'float64'},
    '27':{'short_name': 'empty', 'long_name':'empty', 'units':'NA', 'type':'NA'},
    '0':{'short_name': 'alti', 'long_name':'Altitude', 'units':'m', 'type':'float64'},
    '71':{'short_name': 'cstar', 'long_name':'CStar Transmissometer', 'units': 'ug/l', 'type':'float64'},
    '61':{'short_name': 'u_def', 'long_name':'user defined', 'units':'V', 'type':'float64'},
    '1000':{'short_name': 'sal', 'long_name':'Salinity (C1 T1)', 'units':'PSU', 'type':'float64'}
}

#because it's easier to not think on long flights and instead write mindless code.
scan_header_lookup = {
    'index':{'long_name': 'Scan Index Number', 'type':'int32'},
    't1_C':{'long_name': 'Temperature 1', 'type':'float64'},
    't2_C':{'long_name': 'Temperature 2', 'type':'float64'},
    'p_dbar':{'long_name': 'Pressure', 'type':'float64'},
    'c1_mS/cm':{'long_name': 'Conductivity 1', 'type':'float64'},
    'c2_mS/cm':{'long_name': 'Conductivity 2', 'type':'float64'},
    'sal_PSU':{'long_name': 'Salinity (T1C1)', 'type':'float64'},
    'o1_ml/l':{'long_name': 'Oxygen 1', 'type':'float64'},
    'empty1_NA':{'long_name': 'Empty', 'type':'float64'},
    'empty2_NA':{'long_name': 'Empty', 'type':'float64'},
    'u_def1_V':{'long_name': 'User Defined 1', 'type':'float64'},
    'u_def2_V':{'long_name': 'User Defined 2', 'type':'float64'},
    'fluoro_ug/l':{'long_name': 'Fluorometer', 'type':'float64'},
    'alti_m':{'long_name': 'Altimeter', 'type':'float64'},
    'cstar_ug/l'{'long_name': 'CStar Transmissometer', 'type':'float64'},
    'lat_ddeg':{'long_name':'Latitude Decimal Degrees', 'type':'float64'},
    'lon_ddeg': {'long_name':'Longitude Decimal Degrees', 'type':'float64'},
    'new_fix_bool': {'long_name':'New Fix?', 'type':'bool_'},
    'nmea_datetime': {'long_name':'NMEA Date/Time', 'type':'datetime64'},
    'scan_datetime': {'long_name':'Scan Date/Time', 'type':'datetime64'},
    'btl_fire_bool': {'long_name':'Bottle Fire', 'type':'bool_'}
}
