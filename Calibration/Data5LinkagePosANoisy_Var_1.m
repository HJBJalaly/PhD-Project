%%%%%%%%%%%%%%%
% Real Postition of EE for given coordinates (without noise)
% L0=400-3 (mm)
% L1=300+6 (mm)
% L1=300+2 (mm)
% L2=550-5 (mm)
% L3=550+5 (mm)
% Dq1=-4 (deg)
% Dq2=-6 (deg)
%
% Xm=Xr+ 1*randn---> var=1;
%%
Xm=[ -237.6255
     -189.5374
     -142.9715
      -92.8555
      -46.3103
     -339.7014
     -290.2176
     -240.8519
     -188.3762
     -141.6008
      -87.9165
      -41.4594
        9.3838
       54.9215
     -394.1565
     -341.5736
     -293.2728
     -240.5758
     -190.5456
     -141.8306
      -88.0530
      -36.0240
       14.0827
       60.2627
      109.8552
     -443.0413
     -394.8365
     -345.2987
     -293.1684
     -242.4954
     -191.4890
     -142.2910
      -85.9786
      -36.9505
       14.5176
       65.0236
      112.6022
      160.0995
     -445.5511
     -397.3243
     -346.0557
     -297.1167
     -244.1836
     -192.5745
     -142.2882
      -89.2087
      -34.6151
       13.9723
       67.9587
      117.5401
      163.2527
     -452.3227
     -399.0795
     -347.0180
     -295.1122
     -245.0868
     -193.8287
     -143.4395
      -87.0127
      -33.5990
       15.5963
       67.7876
      119.5111
      166.5741
     -452.3112
     -400.6802
     -350.7746
     -298.4387
     -249.5495
     -193.8704
     -144.3011
      -87.7601
      -34.5076
       16.5399
       69.3449
      119.6873
      170.0386
     -505.1562
     -454.7424
     -401.3257
     -351.6018
     -301.5464
     -248.6730
     -196.5904
     -144.2149
      -88.8524
      -32.2944
       18.5629
       73.1558
      120.4770
      170.6059
      221.1203
     -504.4654
     -455.8194
     -405.1217
     -355.4053
     -303.8605
     -249.0719
     -197.6089
     -143.4492
      -88.2821
      -33.7429
       20.8681
       71.3485
      121.7477
      172.3305
      223.3455
     -509.1486
     -458.2092
     -406.9577
     -356.3920
     -304.9920
     -250.6783
     -127.4901
     -146.8802
      -88.3916
      -34.7871
       19.7462
       73.1469
      124.4250
      173.3068
      224.8887
     -509.3309
     -461.6304
     -410.6395
     -357.2221
     -307.2771
     -252.5021
     -199.8801
     -146.2830
      -88.8772
      -31.9906
       19.4748
       72.9189
      126.5753
      176.1073
      225.0525
     -509.1387
     -460.5068
     -409.7437
     -359.6459
     -306.1974
     -254.3652
     -202.0005
     -145.4305
      -89.5416
      -32.4226
       21.0684
       74.5144
      128.7560
      175.7170
      227.6260
     -459.5841
     -408.0948
     -199.7276
     -146.3648
      -88.3242
      124.3659
      176.9919];
  %%
Ym=[  871.5308
      871.1903
      871.2412
      867.1704
      867.0083
      822.3969
      823.1208
      825.2811
      825.5427
      822.9515
      824.0908
      822.7310
      820.6846
      818.9572
      774.7129
      776.7291
      776.7066
      774.1342
      777.3262
      777.5225
      775.6734
      776.8678
      774.7507
      773.2321
      772.4244
      724.4025
      724.7828
      727.0129
      727.9662
      728.7765
      730.3094
      732.3557
      729.2140
      729.8945
      726.6383
      726.3186
      723.7332
      721.7449
      675.4827
      676.6460
      676.8528
      677.9912
      679.0478
      677.5920
      680.4670
      683.7968
      679.7159
      675.7449
      677.1600
      677.3374
      672.8602
      624.8286
      626.5770
      628.2879
      629.6637
      633.0326
      630.7314
      633.5069
      634.1968
      631.1122
      632.2461
      630.9510
      627.9332
      625.6226
      575.3471
      576.3535
      576.6264
      582.5724
      582.2093
      584.7001
      582.4282
      584.8098
      584.0975
      583.7599
      579.0171
      578.8382
      574.8902
      523.3238
      523.9680
      525.8891
      530.3365
      531.2731
      532.5659
      533.7558
      534.6840
      536.0967
      534.8143
      533.4298
      531.9171
      530.1124
      527.9324
      525.2595
      471.1418
      471.6923
      476.6783
      477.2710
      480.2523
      483.7267
      485.9155
      487.0430
      486.8079
      486.5127
      484.4329
      481.9999
      479.2938
      477.0606
      475.4276
      421.0212
      424.0958
      423.5166
      426.1500
      429.9489
      430.7376
      103.6162
      438.0523
      435.5669
      435.0004
      436.5152
      431.7159
      428.0550
      428.4553
      424.1644
      366.5484
      369.4691
      371.3783
      375.4614
      377.7435
      380.3360
      381.5583
      385.4807
      386.4200
      384.4961
      385.6585
      382.0638
      380.5536
      377.0444
      376.2289
      313.2044
      315.2842
      316.1938
      319.3261
      323.9977
      326.7265
      331.7109
      335.3976
      336.3351
      336.9121
      334.3358
      329.9866
      324.5734
      326.2236
      323.1930
      257.0895
      261.9089
      278.6291
      284.9636
      286.0844
      272.2595
      269.9205];