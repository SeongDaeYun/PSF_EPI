# PSF_EPI
This script simulates point-spread-function (PSF) of echo-planar-imaging (EPI) along one phase encoding direction and calculates the corresponding full-width-half-maximum (FWHM). 


## Contents
[1. How to use](#How-to-use) <br>
[2. Example](#Example)


### How to use
```
Function: [FWHM, h0] =
          PSF_EPI_PE( sz_y, T2, TE, espc, pi_f, pf_f, is_T2d, is_PF1 )
```

`Input : ( sz_y, T2, TE, espc, pi_f, pf_f, is_T2d, is_PF1 )`

| Input | Default value | Description |
| ------ | ------ | ------ |
| sz_y   | N/A | Phase encoding size |
| T2     | N/A | T2 or T2* time |
| TE     | N/A | Echo time (ms) |
| espc   | N/A | Echo spacing time (ms) |
| pi_f   | N/A | Parallel imaging acceleration factor, e.g.) 1, 2, 3, ... |
| pf_f   | N/A | Partial Fourier imaging factor ( >= 0.5 ), e.g.) 5/8, 6/8, ... |
| is_T2d | true | Option to enable/disable T2 decay  |
| is_PF1 | true | Option to simulate PF trajectory, i.e. true: one-sided PF, false: two-sided PF |

`Output : [ FWHM, h0 ]`

| Output | Description |
| ------ | ------ |
| FWHM   | calculated FWHM (in pixel) |
| h0     | figure handle |


### Example
```
FWHM = PSF_EPI_PE( 192, 33.2, 35, 1.2, 3, 6/8 ) ;
FWHM =

    1.7421
```


![Figure](https://github.com/SeongDaeYun/PSF_EPI/blob/main/Figure/Fig1.jpg)


