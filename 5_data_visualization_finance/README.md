# Finance_visualization
Developing some Dash app visualization for simple finance calculations such as Interest, SIP, EMI.

# 

[Tool snapshot](https://github.com/vmos1/Code_highlights/blob/main/5_data_visualization_finance/images/snapshot_dashapp.png)

## Building the conda environment 
Required packages : 
- dash
- pandas
- numpy
  
## Formulae
If `P` is the principal, `r` is the rate of interest per term (in percent) and `n` is the number of terms, 
### Simple Interest 
$$ A = P ( 1 + r \times n ) $$

### Compound Interest 
$$ A = P {\left(1 + r \right)}^n $$

### SIP (Systematic Investment Plan) 
Investing amount $ x_0 $ every term,

$$ x_n = x_0 \left( \frac{R^{n+1}-1}{R-1}\right) $$ 

where  R = 1 + r  

### EMI (Estimated monthly/ term payment) 
If `X` is the payment amount per term and $x_0$ is the initial loan amount,

#### Amount per term 
$$ X = x_0 \left( \frac{R-1}{1-R^{-n}} \right) $$ 

### Number of terms 
$$ n =  - \frac{\log{ \left[1-\frac{x_0(R-1)}{X} \right] }}{\log{R}} $$
