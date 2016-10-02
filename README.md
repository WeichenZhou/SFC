# SFC
It is a new method named Standardized Fold Change (SFC) to detect differential expressions for cases and controls.

# Usage
```
SFC -b bin_size -i input_file -control number_of_controls -case number_of_cases
```
-b: The bin/neighborhood size of the method (>=2). The default value is 1000.

-i: One tab delimited input file is reqiured. No missing values are allowed.

-control: The number of controls for SFC (>=1). 

-case: The number of cases for SFC (>=1).

Output: SCF_output.txt contains values of SFC statistic for each probe. If all cases and controls have the same expression value, value of SFC should be 'nan'.

