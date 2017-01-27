%%%%%%%%%%%%%%%
% Real Postition of EE for given coordinates (without noise)
% L0=400+5 (mm)
% L1=300-3 (mm)
% L1=300-4 (mm)
% L2=550+10 (mm)
% L3=550+8 (mm)
% Dq1=-2 (deg)
% Dq2=1 (deg)
%
% Xm=Xr+ 1*randn---> var=1;
%%
Xm=[ -274.8338
     -225.0620
     -181.4704
     -133.2203
      -85.3035
     -375.4245
     -328.1229
     -281.3493
     -230.5737
     -183.8760
     -133.4974
      -84.0390
      -35.5323
       11.8858
     -429.2326
     -383.2663
     -334.8955
     -281.2076
     -233.0263
     -180.1469
     -128.3503
      -79.4773
      -29.1686
       16.8882
       67.5069
     -484.1410
     -436.1858
     -385.4164
     -335.7567
     -284.9382
     -233.5263
     -180.9860
     -126.4312
      -77.7691
      -26.5985
       27.5417
       75.6628
      122.4938
     -486.7443
     -440.1534
     -390.5578
     -338.9101
     -285.8937
     -232.9072
     -183.5818
     -127.3398
      -71.2073
      -21.7490
       29.7045
       81.1447
      128.0136
     -495.2088
     -442.8528
     -390.9784
     -340.1483
     -285.2596
     -232.1844
     -178.3197
     -123.0528
      -70.2829
      -17.5130
       35.2128
       88.0421
      137.1500
     -497.9677
     -447.7776
     -394.8663
     -342.2909
     -288.7683
     -232.0574
     -179.8644
     -120.3092
      -64.7821
      -10.5868
       40.8266
       92.7764
      143.2032
     -552.4230
     -501.5893
     -448.0048
     -398.3839
     -343.1079
     -285.5601
     -230.9626
     -177.6617
     -118.4811
      -59.4562
       -6.6987
       47.4733
       98.9133
      148.9447
      202.0349
     -553.4714
     -503.9711
     -452.6508
     -401.9257
     -344.7302
     -288.0173
     -232.6958
     -173.2288
     -115.2325
      -56.4434
        0.2832
       54.3864
      106.7438
      158.1463
      206.4216
     -558.5178
     -509.1994
     -455.0267
     -401.6044
     -346.8187
     -289.8299
     -116.1210
     -173.9730
     -111.0532
      -51.8795
        5.9255
       63.1349
      112.2430
      163.9480
      213.3662
     -562.5971
     -511.9383
     -461.1873
     -404.4524
     -350.1081
     -291.4740
     -229.8951
     -170.2135
     -105.6880
      -45.2803
       12.8873
       67.1201
      120.4876
      170.6752
      218.8412
     -565.4208
     -513.6252
     -461.1618
     -409.2238
     -353.5986
     -294.3973
     -230.0871
     -167.5637
     -100.5791
      -39.7034
       19.1137
       72.8269
      128.4288
      176.6944
      225.6235
     -516.0853
     -463.6903
     -230.4879
     -165.2674
      -94.4798
      128.9540
      182.5354];
  %%
Ym=[  881.4626
      884.7182
      883.6670
      886.4782
      887.0965
      829.5911
      836.0521
      836.4444
      841.9695
      842.4588
      846.0453
      847.0179
      847.7496
      848.0850
      783.7566
      788.8340
      791.6957
      795.0084
      798.7073
      800.9478
      803.6605
      804.9417
      803.3422
      805.3777
      805.7712
      733.9465
      738.5430
      742.8170
      747.9671
      751.9679
      754.7619
      759.4954
      758.2024
      758.8877
      760.2667
      761.9265
      760.4969
      760.7568
      687.2506
      694.5875
      695.5469
      703.3112
      706.3895
      709.7353
      712.7301
      715.8860
      716.4231
      715.4822
      711.9842
      715.3034
      716.2942
      641.8557
      648.4777
      653.8452
      657.6195
      662.9554
      668.0512
      668.9578
      671.0228
      672.4489
      672.6440
      670.3783
      669.9463
      667.8711
      595.8670
      601.4085
      606.3440
      613.3757
      618.2950
      621.8233
      622.0148
      626.8278
      624.2591
      626.1430
      625.5181
      621.7712
      621.6282
      540.4751
      546.0575
      554.1816
      560.9588
      567.4978
      572.1602
      575.7609
      579.6758
      580.2168
      580.0342
      578.0563
      578.4899
      575.7126
      577.1573
      570.9500
      492.7771
      497.3460
      506.2145
      511.9629
      521.1958
      528.6973
      531.5316
      536.2004
      534.5157
      532.4956
      532.1913
      529.8568
      527.2956
      522.3088
      522.9136
      439.3429
      449.2753
      457.9406
      464.8331
      473.2336
      480.8595
      182.8199
      489.2431
      489.0630
      486.3531
      486.4006
      480.6970
      477.9874
      476.6532
      472.5651
      388.4600
      395.7323
      407.6787
      415.3218
      425.5845
      432.8209
      438.3306
      441.3238
      441.8151
      439.7060
      436.9702
      431.2830
      426.1550
      425.8733
      420.2438
      334.3603
      342.0299
      353.0768
      364.3285
      376.3089
      388.3019
      392.7430
      397.2850
      396.8934
      392.2829
      384.1523
      376.9220
      372.1630
      369.2853
      368.4327
      282.1811
      293.8535
      347.8571
      352.6831
      349.0454
      310.9738
      308.3989];