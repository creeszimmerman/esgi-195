# esgi-195
GFL GFM toy models


'GFM'
2 GFMs
- want to test if the current limiter coming into play leads to instability
- this is a simpler model than the GFLs, as it is just P-w droop control plus some impedances in the circuit

'GFL'
2 GFLs
- this has PI PLL control and PI current control, plus some impedances in a circuit
- can make the current spike higher by reducing resistance of the circuit
- this qualitatively matches the behaviour seen in Fig 17 in the paper
- you can also play with the controller parameters

'GFL_GFM'
GFL + GFM
- this isn't quite perfect yet as it isn't quite equilibrated initially. I was also trying replicate Fig 18 in the paper but not got there yet
