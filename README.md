# EV
This repository includes scripts for data analysis in the EV paper

## Multi-Trait genome-wide association analysis 
This part uses the MulitABEL packages in R, see <https://github.com/xiashen/MultiABEL>

## run PASCAL using command line
Pascal user mannual and relative paper can be found <https://www2.unil.ch/cbg/index.php?title=Pascal>, we used a different reference genome from TWINSUK
```
    ./Pascal --pval=resources/GIANT/giant_file.txt --customdir=resources/TWINSUK.and.ALSPAC --custom=TWINSUK.and.ALSPAC
```
