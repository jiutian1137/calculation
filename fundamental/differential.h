#pragma once
/*{ "clmagic/calculation/fundamental/differential":{ 
"Mathematician":[ 
  "Brook Taylor", 
  "C.J.F. Ridders" 
] ,

"Reference":{
  "Library": "boost::math::differentiation::differentiation" ,
  "Book": "Numerical Analysis" ,
  "Url":"http://www.holoborodko.com/pavel/about/"
},

"Author": "LongJiangnan",
  
"License": "Please identify Mathematician"
} }*/


#include <cmath>// abs, sqrt, pow, nextafter
#include <numeric>// std::numeric_limits<T>::epsilon, std::numeric_limits<T>::max
namespace calculation {
  /* two point forward difference
  * ...
  */
  template<typename Function, typename Real>
  Real f2_difference(Function f, Real x, Real* error = nullptr) {
    Real eps = std::numeric_limits<Real>::epsilon();
    Real h = static_cast<Real>(sqrt(4 * eps));
	if ( (x + h) - x == 0 ) {
      // similar nextafter(x)-x
      h = x * eps;
    }

    Real y0 = f(x);
    Real yh = f(x + h);
    Real dfx = (yh - y0) / h;
    if (error) {
      // math_error: h/2*ddf(x)
      *error = h / 2 * abs(yh - h * dfx - y0/* + dddf(x) + round_error */);
      // round_error
      *error += (abs(yh) * eps + abs(y0) * eps) / h;
    }
    return dfx;
  }

  /* three point centered difference
  * about 'error', we ignore low order terms in 'formula'
  * about 'h', we remove relation between fxxx(x) and Error
  
  * formula
                 f(x)                fx(x)              fxx(x)              fxxx(x)
  (1) f(x+h) = --------*pow(h,0) + --------*pow(h,1) + --------*pow(h,2) + --------*pow(h,3) +  ...    ,apply 'Taylor theorem'
                fact(0)             fact(1)             fact(2)             fact(3)
                 f(x)                 fx(x)                fxx(x)              fxxx(x)
  (2) f(x-h) = --------*pow(-h,0) + --------*pow(-h,1) + --------*pow(-h,2) + --------*pow(-h,3) +  ...    ,apply 'Taylor theorem'
                fact(0)              fact(1)               fact(2)             fact(3)
       f(x+h) - f(x-h)                       fxx(x)              fxxx(x)                                      fxx(x)              fxxx(x)
  (3) ----------------- =  f(x) + fx(x)*h + --------*pow(h,2) + --------*pow(h,3) +  ...  - f(x) + fx(x)*h - --------*pow(h,2) + --------*pow(h,3) + ...
             2*h                             fact(2)             fact(3)                                      fact(2)             fact(3)
                          -------------------------------------------------------------------------------------------------------------------------------
                                                                                2*h
                                         fxxx(x)
                           2*fx(x)*h + 2*-------*pow(h,3) + ...
                                         fact(3)
                        = --------------------------------------
                                      2*h
                                   fxxx(x)
                        = fx(x) + --------*pow(h,2) + ...
                                   fact(3)
       f(x+h) - f(x-h)     fxxx(x)
      ----------------- - --------*pow(h,2) - ... = fx(x)
             2*h           fact(3)

  * computation formula
                      f(x+h)*(1 +- eps) - f(x-h)*(1 +- eps)
  c3_difference(x) = --------------------------------------
                                       2*h
                      f(x+h) +- f(x+h)*eps - f(x-h) +- f(x-h)*eps
                   = ---------------------------------------------
                                          2*h
                      f(x+h) - f(x-h)      f(x+h)*eps +- f(x-h)*eps
                   = ----------------- +- --------------------------
                            2*h                      2*h

  * error formula
  Error = abs(ErrorResult - NoErrorResult)

  * mathematical error
               f(x+h) - f(x-h)     f(x+h) - f(x-h)     fxxx(x)
  Error = abs(----------------- - ----------------- + --------*pow(h,2) + ...) 
                     2*h                 2*h           fact(3)
                fxxx(x)
        = abs( --------*pow(h,2) + ... )

  * computation error
                f(x+h) - f(x-h)     f(x+h) - f(x-h)      f(x+h)*eps +- f(x-h)*eps
  Error = abs( ----------------- - ----------------- +- -------------------------- )
                      2*h                 2*h                       2*h
            abs(f(x+h))*eps + abs(f(x-h))*eps
        <= -----------------------------------
                           2*h
            max( abs(f(x+h)),abs(f(x-h)) ) * (eps + eps)
        <= ----------------------------------------------
                              2*h
            abs(f(x)) * (eps + eps)
        <= ------------------------- ,apply 'intermediate theorem'
                     2*h
        <= abs(f(x)) * eps / h

  * 'h'
                          fxxx(x)
  (1) TotalError <= abs( --------*pow(h,2) + ... ) + abs(f(x))*eps/h
                          fact(3)
                           fxxx(x)
                 ~<= abs( --------*pow(h,2) ) + abs(f(x))*eps/h  ,ignore low order terms
                           fact(3)
                           pow(h,2)
                 ~<= abs( --------- ) + abs(f(x))*eps/h ,remove relatin between fxxx(x) and Error
                           fact(3)
                 ~<= pow(h,2)/6 + abs(f(x))*eps/h
  (2)
      Error| +            +
           |  +        +
           |   +    +
           |     + 
           +--------------------- h
  (2)            d/dh * TotalError == 0   ,find critial point
      h/3 - abs(f(x))*eps/pow(h,2) == 0
	                           h/3 == abs(f(x))*eps/pow(h,2)  ,+ abs(f(x))*eps/pow(h,2)
                        pow(h,3)/3 == abs(f(x))*eps           ,* pow(h,2)
                          pow(h,3) == abs(f(x))*eps*3
                                 h == pow(abs(f(x))*eps*3, 1/3) ,pow 1/3

  * reference 
  <<Numerical Analysis>> Timothy Sauer
  "https://www.uio.no/studier/emner/matnat/math/MAT-INF1100/h08/kompendiet/diffint.pdf"
  */
  template<typename Function, typename Real>
  Real c3_difference(Function f, Real x, Real* error = nullptr) {
    Real eps = std::numeric_limits<Real>::epsilon();
    Real h = pow(abs(f(x))*eps*3, Real(1)/3);
	if ( (x + h) - x == 0 ) {
      // similar nextafter(x)-x
      h = x * eps;
    }

    Real y_h = f(x + h);
    Real y_m_h = f(x - h);
    if (error) {
      // math_error: pow(h,2)/6 * d3f(x), five_point_centered_difference_third_derivative
      Real y_two_h = f(x + 2 * h);
      Real y_m_two_h = f(x - 2 * h);
      *error = abs((y_two_h - y_m_two_h) - 2 * (y_h - y_m_h)) / (12 * h);
      // round_error: ...
      *error += (abs(y_h) * eps + abs(y_m_h) * eps) / (2 * h);
    }
    return (y_h - y_m_h) / (2 * h);
  }

  /* five point centered difference
  */
  template<typename Function, typename Real>
  Real c5_difference(Function f, Real x, Real* error = nullptr) {
    Real eps = std::numeric_limits<Real>::epsilon();
    Real h = pow(abs(f(x))*eps*Real(11.25), Real(1)/5);
	if ( (x + h) - x == 0 ) {
      // similar nextafter(x)-x
      h = x * eps;
    }

    Real y_h = f(x + h);
    Real y_m_h = f(x - h);
	Real y_two_h = f(x + 2 * h);
	Real y_m_two_h = f(x - 2 * h);
	if (error) {
      // math_error: pow(h,4)/30 * dddddf(x), dddddfx=seven_point_centered_difference_fifth_derivative
      Real y_three_h = f(x + 3 * h);
      Real y_m_three_h = f(x - 3 * h);
      *error = abs((y_three_h - y_m_three_h) + 5 * (y_h - y_m_h) - 4 * (y_two_h - y_m_two_h)) / (60 * h);
      // round_error
      *error += (8 * (abs(y_h) + abs(y_m_h)) + (abs(y_two_h) + abs(y_m_two_h))) * eps / (12 * h);
	}
    return (8 * (y_h - y_m_h) - (y_two_h - y_m_two_h)) / (12 * h);
  }

  /* seven point centered difference
  */
  template<typename Function, typename Real>
  Real c7_difference(Function f, Real x, Real* error = nullptr) {
    Real eps = std::numeric_limits<Real>::epsilon();
    Real h = pow(abs(f(x))*eps*385/9, Real(1.0)/7);
	if ( (x + h) - x == 0 ) {
      // similar nextafter(x)-x
      h = x * eps;
    }

    Real y_h = f(x + h);
	Real y_m_h = f(x - h);
	Real y_two_h = f(x + 2 * h);
	Real y_m_two_h = f(x - 2 * h);
	Real y_three_h = f(x + 3 * h);
	Real y_m_three_h = f(x - 3 * h);
	if (error) {
      // math_error: pow(h,6)/140 * d7fx, nine_point_centered_difference_seventh_derivative
      Real y_four_h = f(x + 4 * h);
	  Real y_m_four_h = f(x - 4 * h);
	  *error = abs((y_four_h - y_m_four_h) - 14 * (y_h - y_m_h) + 14 * (y_two_h - y_m_two_h) - 6 * (y_three_h - y_m_three_h)) / (280 * h);
	  // round_error: ...
	  *error += (abs(y_three_h) + abs(y_m_three_h) + 9 * (abs(y_two_h) + abs(y_m_two_h)) + 45 * (abs(y_h) + abs(y_m_h))) * eps / (60 * h);
	}
	return ((y_three_h - y_m_three_h) - 9 * (y_two_h - y_m_two_h) + 45 * (y_h - y_m_h)) / (60 * h);
  }

  /* nine point centered difference
    Notation:{
			f(x) = f(x),
			fx(x) = dfdx(x),
			fxx(x) = ddfdxx(x) = fx2(x),
			fxxx(x) = dddfdxxx(x) = fx3(x),
			fxxxx(x) = ddddfdxxxx(x) = fx4(x),
			fx...x(x) = d..dfdx..x(x) = fxN(x)
		},
		Taylor-theorem:{
			f(x+h) = pow(h,0) * (1/fact(0)) * f(x)
				   + pow(h,1) * (1/fact(1)) * fx(x)
				   + pow(h,2) * (1/fact(2)) * fxx(x)
				   + ... 
				   + pow(h,N-1) * (1/fact(N-1)) * df(N-1x)(x)
				   + integrate[ pow(x+h-u,N-1) * (1/fact(N-1)) * dN-1f(u) * du ]
		},
		Step1:{ "Use Taylor theorem, and ignore unwanted terms",
			f(x+h) =  pow(h,0) * (1/1) * f(x)
					+ pow(h,1) * (1/1) * fx(x)
					+ pow(h,2) * (1/2) * fxx(x)
					+ pow(h,3) * (1/6) * fxxx(x)
					+ pow(h,4) * (1/24) * fxxxx(x)
					+ pow(h,5) * (1/120) * fxxxxx(x)
					+ pow(h,6) * (1/720) * fxxxxxx(x)
					+ pow(h,7) * (1/5'040) * fxxxxxxx(x)
					+ pow(h,8) * (1/40'320) * fxxxxxxxx(x)
					+ pow(h,9) * (1/362'880) * fxxxxxxxxx(x)
					+ pow(h,10) * (1/3'628'800) * fxxxxxxxxxx(x)
					+ pow(h,11) * (1/39'916'800) * fxxxxxxxxxxx(x)
					+ ...ignore;
			f(x-h) = ...;

			f(x+2*h) =  pow(2*h,0) * (1/1) * f(x)
					  + pow(h,1) * (2/1) * fx(x)
					  + pow(h,2) * (4/2) * fxx(x)
					  + pow(h,3) * (8/6) * fxxx(x)
					  + pow(h,4) * (16/24) * fxxxx(x)
					  + pow(h,5) * (32/120) * fxxxxx(x)
					  + pow(h,6) * (64/720) * fxxxxxx(x) 
					  + pow(h,7) * (128/5'040) * fxxxxxxx(x)
					  + pow(h,8) * (256/40'320) * fxxxxxxxx(x)
					  + pow(h,9) * (512/362'880) * fxxxxxxxxx(x)
					  + pow(h,10) * (1'024/3'628'800) * fxxxxxxxxxx(x)
					  + pow(h,11) * (2'048/39'916'800) * fxxxxxxxxxxx(x)
					  + ...ignore;
			f(x-2*h)) = ...;

			f(x+3*h) =  pow(3*h,0) * (1/1) * f(x)
					  + pow(h,1) * (3/1) * fx(x)
					  + pow(h,2) * (9/2) * fxx(x)
					  + pow(h,3) * (27/6) * fxxx(x)
					  + pow(h,4) * (81/24) * fxxxx(x)
					  + pow(h,5) * (243/120) * fxxxxx(x)
					  + pow(h,6) * (729/720) * fxxxxxx(x)
					  + pow(h,7) * (2'187/5'040) * fxxxxxxx(x)
					  + pow(h,8) * (8'748/40'320) * fxxxxxxxx(x)
					  + pow(h,9) * (19'683/362'880) * fxxxxxxxxx(x)
					  + pow(h,10) * (59'049/3'628'800) * fxxxxxxxxxx(x)
					  + pow(h,11) * (177'147/39'916'800) * fxxxxxxxxxxx(x)
					  + ...ignore;
			f(x-3*h) = ...;

			f(x+4*h) =  pow(4*h,0) * (1/1) * f(x)
					  + pow(h,1) * (4/1) * fx(x)
					  + pow(h,2) * (16/2) * fxx(x)
					  + pow(h,3) * (64/6) * fxxx(x)
					  + pow(h,4) * (256/24) * fxxxx(x)
					  + pow(h,5) * (1'024/120) * fxxxxx(x)
					  + pow(h,6) * (4'096/720) * fxxxxxx(x)
					  + pow(h,7) * (16'384/5'040) * fxxxxxxx(x)
					  + pow(h,8) * (65'536/40'320) * fxxxxxxxx(x)
					  + pow(h,9) * (262'144/362'880) * fxxxxxxxxx(x)
					  + pow(h,10) * (1'048'576/3'628'800) * fxxxxxxxxxx(x)
					  + pow(h,11) * (4'194'304/39'916'800) * fxxxxxxxxxxx(x)
					  + ...ignore;
			f(x-4*h) = ...;
		},
		Step2:{ "Get dfdx(x) = ?, use elimination and sorting",
			f(x+h)-f(x-h) =  
				  pow(h,1) * (2/1) * fx(x)
				+ pow(h,3) * (1/3) * fx3(x)
				+ pow(h,5) * (1/60) * fx5(x)
				+ pow(h,7) * (1/2'520) * fx7(x)
				+ pow(h,9) * (1/181'440) * fx9(x)
				+ pow(h,11) * (1/19'958'400) * fx11(x)
				+ ...ignore;

			f(x+2*h)-f(x-2*h) = 
				  pow(h,1) * (4/1) * fx(x)
				+ pow(h,3) * (8/3) * fx3(x)
				+ pow(h,5) * (32/60) * fx5(x)
				+ pow(h,7) * (128/2'520) * fx7(x)
				+ pow(h,9) * (512/181'440) * fx9(x)
				+ pow(h,11) * (2'048/19'958'400) * fx11(x)
				+ ...ignore;

			f(x+3*h)-f(x-3*h) = 
				  pow(h,1) * (6/1) * fx(x)
				+ pow(h,3) * (27/3) * fx3(x)
				+ pow(h,5) * (243/60) * fx5(x)
				+ pow(h,7) * (2'187/2'520) * fx7(x)
				+ pow(h,9) * (19'683/181'440) * fx9(x)
				+ pow(h,11) * (177'147/19'958'400) * fx11(x)
				+ ...ignore;

			f(x+4*h)-f(x-4*h) = 
				  pow(h,1) * (8/1) * fx(x)
				+ pow(h,3) * (64/3) * fx3(x)
				+ pow(h,5) * (1'024/60) * fx5(x)
				+ pow(h,7) * (16'384/2'520) * fx7(x)
				+ pow(h,9) * (262'144/181'440) * fx9(x)
				+ pow(h,11) * (4'194'304/19'958'400) * fx11(x)
				+ ...ignore;

			672 * [ f(x+h) - f(x-h) ] =
				  pow(h,1) * (672*2) * fx(x)
				+ pow(h,3) * (672/3) * fx3(x)
				+ pow(h,5) * (672/60) * fx5(x)
				+ pow(h,7) * (672/2'520) * fx7(x)
				+ pow(h,9) * (672/181'440) * fx9(x)
				+ pow(h,11) * (672/19'958'400) * fx11(x)
				+ ...ignore;

			168 * [ f(x+2*h) - f(x-2*h) ] =
					pow(h,1) * (168*4) * fx(x)
				+ pow(h,3) * (168*8 / 3) * fx3(x)
				+ pow(h,5) * (168*32 / 60) * fx5(x)
				+ pow(h,7) * (168*128 / 2'520) * fx7(x)
				+ pow(h,9) * (168*512 / 181'440) * fx9(x)
				+ pow(h,11) * (168*2'048 / 19'958'400) * fx11(x)
				+ ...ignore;

			32 * [ f(x+3*h) - f(x-3*h) ] =
					pow(h,1) * (32*6) * fx(x)
				+ pow(h,3) * (32*27 / 3) * fx3(x)
				+ pow(h,5) * (32*243 / 60) * fx5(x)
				+ pow(h,7) * (32*2'187 / 2'520) * fx7(x)
				+ pow(h,9) * (32*19'683 / 181'440) * fx9(x)
				+ pow(h,11) * (32*177'147 / 19'958'400) * fx11(x)
				+ ...ignore;

			3 * [ f(x+4*h) - f(x-4*h) ] =
					pow(h,1) * (3*8) * fx(x)
				+ pow(h,3) * (3*64 / 3) * fx3(x)
				+ pow(h,5) * (3*1'024 / 60) * fx5(x)
				+ pow(h,7) * (3*16'384 / 2'520) * fx7(x)
				+ pow(h,9) * (3*262'144 / 181'440) * fx9(x)
				+ pow(h,11) * (3*4'194'304 / 19'958'400) * fx11(x)
				+ ...ignore;

			672*[f(x+h)-f(x-h)] - 168*[f(x+2*h)-f(x-2*h)] + 32*[f(x+3*h)-f(x-3*h)] - 3*[f(x+4*h)-f(x-4*h)]
				=    pow(h,1) * (672*2 - 168*4 + 32*6 - 3*8) * fx(x)
				+ pow(h,3) * (672/3 - 168*8/3 + 32*27/3 - 3*64/3) * fx3(x)
				+ pow(h,5) * (672/60 - 168*32/60 + 32*243/60 - 3*1'024/60) * fx5(x)
				+ pow(h,7) * (672/2'520 - 168*128/2'520 + 32*2'187/2'520 - 3*16'384/2'520) * fx7(x)
				+ pow(h,9) * (672/181'440 - 168*512/181'440 + 32*19'683/181'440 - 3*262'144/181'440) * fx9(x)
				+ pow(h,11) * (672/19'958'400 - 168*2'048/19'958'400 + 32*177'147/19'958'400 - 3*4'194'304/19'958'400) * fx11(x)
				+ ...ignore
				=    pow(h,1) * 840 * fx(x)
				+ pow(h,3) * 0 * fx3(x)
				+ pow(h,5) * 0 * fx5(x)
				+ pow(h,7) * 0 * fx7(x)
				+ pow(h,9) * (-4/3) * fx9(x)
				+ ...ignore
				= 840*h*fx(x) - (4/3)*pow(h,9)*fx9(x) + ...ignore;

			fx(x) * 840*h = { 672*[f(x+h)-f(x-h)] - 168*[f(x+2*h)-f(x-2*h)] + 32*[f(x+3*h)-f(x-3*h)] - 3*[f(x+4*h)-f(x-4*h)] }
							+ (4/3) * pow(h,9) * fx9(x)
							+ ...ignore
			fx(x)		  = { 672*[f(x+h)-f(x-h)] - 168*[f(x+2*h)-f(x-2*h)] + 32*[f(x+3*h)-f(x-3*h)] - 3*[f(x+4*h)-f(x-4*h)] } / (840*h)
				            + (4/3) * pow(h,9) * fx9(x) / (840*h) 
							+ ...igonore / (840*h) 
			
			dfdx(x)		  = { 672*[f(x+h)-f(x-h)] - 168*[f(x+2*h)-f(x-2*h)] + 32*[f(x+3*h)-f(x-3*h)] - 3*[f(x+4*h)-f(x-4*h)] } / (840*h)
							+ (1/630) * pow(h,8) * fx9(x)
							+ ...igonore
		}
	} */
  template<typename Function, typename Real>
  Real c9_difference(Function f, Real x, Real* error = nullptr) {
    Real eps = std::numeric_limits<Real>::epsilon();
    Real h = pow(abs(f(x))*eps*2625/16, Real(1)/9);
	if ( (x + h) - x == 0 ) {
      // similar nextafter(x)-x
      h = x * eps;
    }

    Real y_h = f(x + h);
    Real y_m_h = f(x - h);
	Real y_two_h = f(x + 2 * h);
	Real y_m_two_h = f(x - 2 * h);
	Real y_three_h = f(x + 3 * h);
	Real y_m_three_h = f(x - 3 * h);
	Real y_four_h = f(x + 4 * h);
	Real y_m_four_h = f(x - 4 * h);
    if (error) {
      // round_error: ...
      *error = (
        672 * (abs(y_h) + abs(y_m_h)) +
        168 * (abs(y_two_h) + abs(y_m_two_h)) +
        32 * (abs(y_three_h) + abs(y_m_three_h)) +
        3 * (abs(y_four_h) + abs(y_m_four_h))) * eps / (840 * h);
      // math_error: pow(h,8)*(4/3)*d9f(x) 
    }
    return (672 * (y_h - y_m_h)
            - 168 * (y_two_h - y_m_two_h)
            + 32 * (y_three_h - y_m_three_h)
            - 3 * (y_four_h - y_m_four_h)) / (840 * h);
  }

  /* derivative of Function('f') respect to 'x'
  
  * differentiate( sin, x )
   = cos(x)

  * differentiate( cos, x )
   = -sin(x)

  ...

  */
  template<typename Function, typename Real>
  Real differentiate(Function f, Real x, Real* error = nullptr) {
    return c5_difference(f, x, error);
  }

  /* del(v)
  */
  template<typename Function, typename Vector3>
  Vector3 gradient(Function f, Vector3 v) {
    using Real = std::remove_cvref_t< decltype( v[0] ) >;
    Real x = v[0];
    Real y = v[1];
    Real z = v[2];
    auto fdx = [f,y,z](Real x){ return f(Vector3{x,y,z}); };
    auto fdy = [f,x,z](Real y){ return f(Vector3{x,y,z}); };
    auto fdz = [f,x,y](Real z){ return f(Vector3{x,y,z}); };
    Real dfdx = differentiate(fdx, x);
    Real dfdy = differentiate(fdy, y);
    Real dfdz = differentiate(fdz, z);
    return Vector3{ dfdx, dfdy, dfdz };
  }

  /* dot(del,v)
  */
  template<typename Function, typename Vector3>
  auto divergence(Function f, Vector3 v) {
    using Real = std::remove_cvref_t< decltype( v[0] ) >;
    Real x = v[0];
    Real y = v[1];
    Real z = v[2];
    Real dfdx_x = differentiate( [f,y,z](Real x){ return f(Vector3{x,y,z})[0]; }, x );
    Real dfdy_y = differentiate( [f,x,z](Real y){ return f(Vector3{x,y,z})[1]; }, y );
    Real dfdz_z = differentiate( [f,x,y](Real z){ return f(Vector3{x,y,z})[2]; }, z );
    return dfdx_x + dfdy_y + dfdz_z;
  }

  /* cross(del,v) 
  */
  template<typename Function, typename Vector3>
  Vector3 curl(Function f, Vector3 v) {
    using Real = std::remove_cvref_t< decltype(v[0]) >;
    Real x = v[0];
    Real y = v[1];
    Real z = v[2];
    Real dfdy_z = differentiate( [f,x,z](Real y){return f(Vector3{x,y,z})[2];}, y );
    Real dfdz_y = differentiate( [f,x,y](Real z){return f(Vector3{x,y,z})[1];}, z );
    Real dfdz_x = differentiate( [f,x,y](Real z){return f(Vector3{x,y,z})[0];}, z );
    Real dfdx_z = differentiate( [f,y,z](Real x){return f(Vector3{x,y,z})[2];}, x );
    Real dfdx_y = differentiate( [f,y,z](Real x){return f(Vector3{x,y,z})[1];}, x );
    Real dfdy_x = differentiate( [f,x,z](Real y){return f(Vector3{x,y,z})[0];}, y );
    return Vector3{ dfdy_z - dfdz_y, dfdz_x - dfdx_z, dfdx_y - dfdy_x };
  }


  // undeterminant ...

  /* three point centered difference, second order derivative
  */
  template<typename Function, typename Real>
  Real c3_difference2nd(Function f, Real x, Real* error = nullptr) {
    Real y = f(x);

    Real eps = std::numeric_limits<Real>::epsilon();
    Real h = static_cast<Real>(pow(48 * abs(y) * eps, 1.0 / 4.0));
    if ( (x + h) - x == 0 ) {
      // similar nextafter(x)-x
      h = x * eps;
    }

    Real y_h = f(x + h);
    Real y_m_h = f(x - h);
    if (error) {
      // math_error: pow(h,2)*(1/12)*d4f(x)
    }
    return (y_h + y_m_h - 2 * y) / pow(h, 2);
  }

  /* five point centered difference, fourth order derivative
  */
  template<typename Function, typename Real>
  Real differentiate4th(Function f, Real x, Real* error = nullptr) {
    Real y = f(x);

    // E = pow(h,2)*(1/6)*d6f(x) + 16*e/pow(h,4)
    // dE = h*(2/6)*d6f(x) - 4*16*e*pow(h,-5) = 0
    // h*(2/6)*d6f(x) = 4*16*e*pow(h,-5)
    // pow(h,6) = 64*(6/2)*e/d6f(x)
    // h approxi pow(192*e, 1.0/6.0)
    Real eps = std::numeric_limits<Real>::epsilon();
    Real h = static_cast<Real>(pow(192 * abs(y) * eps, 1.0 / 6.0));
    if ( (x + h) - x == 0 ) {
      // similar nextafter(x)-x
      h = x * eps;
    }

    Real y_h = f(x + h);
	Real y_m_h = f(x - h);
	Real y_two_h = f(x + 2 * h);
	Real y_m_two_h = f(x - 2 * h);
	if (error) {
      // math_error: pow(h,2)*(1/6)*d6f(x)
	}
	return ((y_two_h + y_m_two_h) - 4 * (y_h + y_m_h) + 6 * y) / pow(h, 4);
  }
}
