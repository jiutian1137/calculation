#pragma once
/*{ "clmagic/calculation/fundamental/integral":{ 
"Article":{
  "Rimemann Sum":{ 
    "Mathematician":"Georg Friedrich Bernhard Riemann" },
    
  "Interpolating Polynomial Integration":{
    "Mathematician":[
      "Isaac Newton",
      "Roger Cotes",
      "Thomas Simpson",
      "George Boole",
      "Abramowitz ...",
      "Stegun ...",
      "Ueberhuber ..."
    ], "Formula":"https://mathworld.wolfram.com/Newton-CotesFormulas.html"
  },
    
  "Gauss Integration":{
    "Mathematician":[
      "Johann Carl Friedrich Gauss"
    ], "Formula":"https://mathworld.wolfram.com/GaussianQuadrature.html"
  },
    
  "MonteCarlo Integration":{

  }
},

"Reference":[
  { "Book":"Numerical Analysis",
    "Author":"Timothy Sauer" },

  { "Book":"Harmonia Mensurarum",
    "Author":"..." },

  { "URL":"https://people.sc.fsu.edu/~sshanbhag/NumericalIntegration.pdf", 
    "Author":"...",
    "Desc":"some example" },

  { "URL":"http://homepage.math.uiowa.edu/~whan/3800.d/S5-2.pdf",
    "Author":"?",
    "Desc":"correct_polynomial_x" }
],

"Author":"LongJiangnan",

"License":"Please Identify Mathmatician"
} }*/


#include <cmath>
#include <numeric>
#include <vector>
namespace calculation {
#ifndef __CLMAGIC_CALULATION_SUM
#define __CLMAGIC_CALULATION_SUM
	template<typename Integer, typename Function>
	auto sum(Integer start, Integer end, Function f) -> decltype(f(start) + f(end)) {
		auto result = f(start);
		for (Integer i = start + 1; i <= end; ++i) {
			result += f(i);
		}

		return std::move(result);
	}
#endif

  /* |   |   |   |   |     |   |
     x   x   x   x   x ... x   |
     |   |   |   |   |     |   | 
  */
  template<typename Real, typename Function>
  auto upper_sum(Real a, Real b, Function f, size_t n) -> decltype(f(a)) {
    Real dx = (b - a) / n;
	
    auto result = static_cast<decltype(f(a))>(0);
    for (size_t i = 0; i != n; ++i) {
      result += f(a + i*dx) * dx;
    }
    return result;
  }

  /* |   |   |   |   |     |   |
     |   x   x   x   x ... x   x
     |   |   |   |   |     |   | 
  */
  template<typename Real, typename Function>
  auto lower_sum(Function f, Real a, Real b, size_t n) -> decltype(f(a)) {
    Real dx = (b - a) / n;
	
    auto result = static_cast<decltype(f(a))>(0);
    for (size_t i = 0; i != n; ++i) {
      result += f(a + (i+1)*dx) * dx;
    }
    return result;
  }

  /* |   |   |   |   |     |   |
     | x | x | x | x | ... | x | 
     |   |   |   |   |     |   | 
  */
  template<typename Real, typename Function>
  auto middle_sum(Real a, Real b, Function f, size_t n) -> decltype(f(a)) {
    Real dx = (b - a) / n;

    auto result = static_cast<decltype(f(a))>(0);
    for (size_t i = 0; i != n; ++i) {
      result += f( a + (i*dx + (i+1)*dx)/2 ) * dx;
    }
    return result;
  }


  /* first order interpolating polynomial composite integration( Trapozoidal )
  * formula: integral( 'a', 'b', (y0 + y1)*(h/2) ) + Error
  * ~= sum( '0', 'n-1', (y0 + y1)*(h/2) )
  * ~= ( f('a') + f('b') + sum(1, ('n'+1)-1, f('a'+dx*i)) * 2 ) * h/2
  * a-----|-----|-----b
  * |--dx-|--h--|
  * 
  * interpolating polynomial: in [x0,x1], f(x) = y0*(x-x1)/(x0-x1) + y1*(x-x0)/(x1-x0) + Error
  * integral( x0,x1, f(x) )
  * = integral( x0,x1, y0*(x-x1)/(x0-x1) + y1*(x-x0)/(x1-x0) + Error )
  * = y0*((x1-x0)/2) + y1*((x1-x0)/2) + Error
  * = trapozoidal_area + Error
  */
  template<typename Real, typename Function>
  auto interp1st_integral(Real a, Real b, Function f, size_t segments = 128) -> decltype(f(a)) {
    Real dx = (b - a) / segments;
    return ( 
        ( f(a) + f(b) )/2 
      + sum(size_t(1),(segments+1)-2,[&](size_t i){ return f(a+dx*i); })/* sum_y0y2_overlap/2 */ ) * dx;
  }

  /* second order interpolating polynomial composite integration( Simpson )
  * formula: integral( 'a', 'b', (y0 + y1*4 + y2)*h/3 ) + Error
  * ~= sum( 0, 'n-1', (y0 + y1*4 + y2)*h/3 )
  * ~= ( f('a') + sum_y1*4 + f('b') + sum_y0y2_overlap ) * h/3
  * a---+---|---+---|---+---b
  * |---dx--|-h-|
  */
  template<typename Real, typename Function>
  auto interp2nd_integral(Real a, Real b, Function f, size_t segments = 96) -> decltype( f(a) ) {
    Real dx = (b - a) / segments;
    Real h = dx / 2;
    return (  
        f(a)
      + sum( size_t(0), segments-1, [&](size_t i){ return f(a + h + dx*i); } ) * 4 /* sum_y1 * 4 */
      + f(b)
      + sum( size_t(1), (segments+1)-2, [&](size_t i){ return f(a + dx*i); } ) * 2 /* sum_y0y2_overlap */ ) * h/3;
  }

  /* third order interpolating polynomial composite integration( Simpson 3/8 )
  * formula: integral( 'a', 'b', (y0 + y1*3 + y2*3 + y3)*h*3/8 ) + Error
  * ~= sum( '0', 'n-1', (y0 + y1*3 + y2*3 + y3)*h*3/8 )
  * ~= ( f('a') + sum_y1*3 + sum_y2*3 + f('b') + sum_y0y3_overlap ) * h*3/8
  * a---+---+---|---+---+---|---+---+---b
  * |-----dx----|-h-|
  */
  template<typename Real, typename Function>
  auto interp3rd_integral(Real a, Real b, Function f, size_t n = 64) -> decltype( f(a) ) {
    Real dx = (b - a) / n;
	Real h = dx / 3;
    return (  
        f(a)
      + sum( size_t(0), n-1, [&](size_t i){ return f(a + h + dx*i); } ) * 3 /* sum_y1 * 3 */
      + sum( size_t(0), n-1, [&](size_t i){ return f(a + h*2 + dx*i); } ) * 3 /* sum_y2 * 3 */
      + f(b)
      + sum( size_t(1), (n+1)-2, [&](size_t i){ return f(a + dx*i); } ) * 2 /* sum_y0y2_overlap */ ) * h*3/8;
  }

  /* fourth order interpolating polynomial composite integration( Boole )
  * formula: integral( 'a', 'b', (y0*7 + y1*32 + y2*12 + y3*32 + y4*7)*h*2/45 ) + Error
  * ~= sum( 0, 'n-1', (y0*7 + y1*32 + y2*12 + y3*32 + y4*7)*h*2/45 )
  * ~= ( f('a')*7 + sum_y1*32 + sum_y2*12 + sum_y3*32 + f('b')*7 + sum_y0y4_overlap*7 ) * h*2/45
  * a---+---+---+---|---+---+---+---|---+---+---+---b
  * |------dx-------|-h-|
  */
  template<typename Real, typename Function>
  auto interp4th_integral(Real a, Real b, Function f, size_t n = 52) -> decltype( f(a) ) {
    Real dx = (b - a) / n;
    Real h = dx / 4;
    return (  
        f(a) * 7
      + sum( size_t(0), n-1, [&](size_t i){ return f(a + h + dx*i); } ) * 32 /* sum_y1 * 32 */
      + sum( size_t(0), n-1, [&](size_t i){ return f(a + h*2 + dx*i); } ) * 12 /* sum_y2 * 12 */
      + sum( size_t(0), n-1, [&](size_t i){ return f(a + h*3 + dx*i); } ) * 32 /* sum_y3 * 32 */
      + f(b) * 7
      + sum( size_t(1), (n+1)-2, [&](size_t i){ return f(a + dx*i); } ) * 2 * 7 /* sum_y0y4_overlap * 7 */ ) * h*2/45;
  }

  /* seventh order interpolating polynomial composite integration( Weddle )
  * formula: integral( 'a', 'b', (y0 + y1*5 + y2 + y3*6 + y4 + y5*5 + y6)*h*3/10 )
  * multipliers are integers, the only divisor is also integers, and very stable
  */
  template<typename Real, typename Function>
  auto weddle_integral(Real a, Real b, Function f, size_t n = 16) {
    Real dx = (b - a) / n;
    Real h = dx / 6;
    return (  
        f(a)
      + sum( size_t(0), n-1, [&](size_t i){ return f(a + h + dx*i); } ) * 5 /* sum_y1 * 5 */
      + sum( size_t(0), n-1, [&](size_t i){ return f(a + h*2 + dx*i); } ) /* sum_y2 */
      + sum( size_t(0), n-1, [&](size_t i){ return f(a + h*3 + dx*i); } ) * 6 /* sum_y3 * 6 */
      + sum( size_t(0), n-1, [&](size_t i){ return f(a + h*4 + dx*i); } ) /* sum_y4 */
      + sum( size_t(0), n-1, [&](size_t i){ return f(a + h*5 + dx*i); } ) * 5 /* sum_y5 * 5 */
      + f(b)
      + sum( size_t(1), (n+1)-2, [&](size_t i){ return f(a + dx*i); } ) * 2 /* sum_y0y6_overlap */ ) * h*3/10;
  }

  /* integral['a' -> 'b']('f'(x), 'dx'|'n')
  
  * integral( u, dv )
   = u*v - integral( v, du )
  
  * integral( pow(x,N), dx )
   = pow(x, N + 1)/(N + 1) ,apply anti-derivative: derivative( pow(x,N) ) = N*pow(x,N-1)
  
  * integral( 1/x, dx )
   = ln(abs(x)) ,apply anti-derivative: derivative( ln(x) ) = 1/x
  
  * integral( sin(x), dx )
   = -cos(x) ,apply anti-derivative: derivative( cos(x) ) = -sin(x)

  * integral( cos(x), dx )
   = sin(x) ,apply anti-derivative: derivative( sin(x) ) = cos(x)

  * integral( sin(x*A), dx )
   = integral( sin(x*A)*-A/-A, dx )
   = integral( sin(x*A)*-A, dx ) / -A
   = cos(x*A) / -A ,apply anti-derivative: derivative( cos(x*A) ) = -sin(x*A)*A
   
  * integral( 1/(x*A + B), dx )
   = ln(abs(x*A + B))/A
  
  * integral( 1/pow(x + A, 2), dx ) = - 1/(x + A)
  
  * integral( pow(x + A, N), dx ) = pow(x + A, N + 1)/(N + 1)

  */
  template<typename Real, typename Function>
  auto integral(Real a, Real b, Function f, size_t n = 16) {
    return weddle_integral(a, b, f, n);
  }

  


  /*gauss quadrature
  * formula: integral( f(x) * dx ) ~= sum( wi * f(xi) ) 
  * domain: [0,1]
  */
  template<typename Function, typename Real>
  auto gauss_integral(Function f, const Real* x_array, const Real* w_array, size_t array_size) {
    auto result = f(x_array[0]) * w_array[0];
    for (size_t i = 1; i != array_size; ++i) {
      result += f(x_array[i]) * w_array[i];
    }

    return result;
  }

  /*gauss quadrature
  * formula: integral( f(x) * dx ) ~= sum( wi * f(xi) ) 
  * domain: [a,b]
  */
  template<typename Function, typename Real>
  auto gauss_integral(Real a, Real b, Function f, const Real* x_array, const Real* w_array, size_t array_size) {
    Real scale = (b - a) / 2;
    Real center = (a + b) / 2;
    auto result = f(x_array[0] * scale + center) * w_array[0];
    for (size_t i = 1; i != array_size; ++i) {
      result += f(x_array[i] * scale + center) * w_array[i];
    }

    return result * scale;
  }

  template<typename Real>
  struct GaussQuadrature {
    std::vector<Real> x_array;
    std::vector<Real> w_array;
    
    /* integral interface */
    template<typename Function>
    auto operator()(Real a, Real b, Function f) const {
      assert( x_array.size() == w_array.size() );
      return gauss_integral(a, b, f, x_array.data(), w_array.data(), x_array.size());
    }
  };
}// namespace calculation