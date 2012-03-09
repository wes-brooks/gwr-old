Tract Broadband Data / Census Demographics Readme

The FCC releases broadband data on a 0-5 scale, with this explanation:

Code	Connections per 1,000 Households
0	0
1	Zero < x <= 200
2	200 < x <= 400
3	400 < x <= 600
4	600 < x <= 800
5	800 < x

The data is released at two standards. One (hhs_MMYY) is the ability to download OR upload data at a rate of 25 KB/second (200 Kbps). The other, identified (btop_MMYY), requires the ability to download data at a rate of 96 KB/second (768 Kbps) and upload files at a rate of 25KB/second (200 Kbps).

The raw data is available for download at the FCC's site here: http://transition.fcc.gov/wcb/iatd/comp.html in the "Census Tract Information Mapped for Internet Access Services faster than 200 kbps in at least one direction" section, although the header is slightly misleading--the files include both the hhs and btop standards. 
