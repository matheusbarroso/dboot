# dboot
dependent bootstrap

This is a package designed to provide plug-in estimates for the optimal block length estimate in the context of the Moving Block Bootstrap. Currently, only the algorithm of Hall, Horowitz and Jign (1995) is working properly. The NPPI algorithm of Lahiri, Furukawa and Lee (2007) is in final development fase. 

## To *install* the package:
Inside **R** console just enter these two line commands:

**library**(devtools) 

**install_github**("matheusbarroso/dboot") 


Check function's help for examples...

References: 

Hall, P., Horowitz, J. L e Jing, B-Y. 1995. On blocking rules for the bootstrap with dependent data. Biometrika. 82, 1995, pp. 561-574.

Lahiri, Soumendra, Furukawa, Kyoji and Lee, Yd. 2007. A nonparametric plug-in rule for selecting optimal block lengths for block bootstrap methods. Statistical Methodology. July, 2007, Vol. 4, 3, pp. 292-321.
 
