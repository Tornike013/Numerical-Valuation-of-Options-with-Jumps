# Numerical-Valuation-of-Options-with-Jumps
Group Assingment in Numerical Valuation of Options with Jumps in the Underlying - Finite Difference Method for Option Pricing under the Merton Jump-Diffusion Model

Introduction:
Various alternative models have been developed in the financial literature to
address issues like stochastic volatility, deterministic local volatility functions,
jump-diffusion models, and L´evy models. Jump-diffusion and L´evy models are
particularly attractive due to their ability to explain stock jump patterns and
realistic pricing of options close to maturity.

The paper focuses on the numerical valuation of European Vanilla options
using the jump-diffusion approach with constant coefficients. It explores models
by Merton and Kou, utilizing analytical formulas for solutions. The research
aims to lay the groundwork for solving more complex financial products like
American options and path-dependent options.

The numerical valuation of jump-diffusion processes involves various techniques such as the ADI finite difference method combined with the fast Fourier
transform and a finite element method. The paper discusses the discretization
process, the solver methodology, and the use of the fast Fourier transform to
enhance computational efficiency. Numerical tests are conducted to verify the
accuracy of the numerical schemes.

This project aims to develop and implement numerical methods for valuing
European Vanilla options under the jump-diffusion framework with constant
coefficients, specifically focusing on the Merton.
While the paper explores analytical solutions, the project emphasizes creating numerical implementations using the techniques described in sections 3 and
4. This can involve:

Reproducing the techniques: The project can replicate the numerical methods presented in section 4, such as the BDF2 finite difference method combined
with the fast Fourier transform (FFT) or a finite element method. Presenting
alternative methods: The project can explore alternative numerical approaches
besides those mentioned in the paper for valuing the options. Section 7 likely
provides reference results obtained using the paper’s methods. The project
should compare the results obtained from the implemented numerical methods
with these reference values to assess their accuracy and efficiency.

In essence, the project aims to bridge the gap between theoretical models
and practical applications by providing functional numerical tools for option
valuation under jump-diffusion models
