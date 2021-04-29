#pragma once
/*{ "clmagic/calculation/formulas/fractional":{
  "Article":{
    "Perlin noise":{ "Author":"Ken Perlin" } ,
    "Simplex noise":{ "Author":["Ashima Arts","Stefan Gustavson"], "Code":"https://github.com/stegu/webgl-noise" } ,
    "Classic noise":{ "Author":"Stefan Gustavson", "Code":"https://github.com/stegu/webgl-noise" } ,
    "Cellular noise":{ "Author":"Stefan Gustavson",  "Code":"https://github.com/stegu/webgl-noise" }
  },
  "Book":{
    "Texturing and Modeling":{ "Author":[], "Date":"2005" }
  },
  "License": "Please Identify Article"
} }*/


#include <type_traits>
#include "../fundamental/array.h"
namespace calculation {
	template<typename Real, typename G> inline 
	Real noise1(Real x, G lattice) {
		// setup
		const Real xi = floor(x);
		const Real xf = x - xi;
		// value = sample lattice
		const Real w0 = lattice(xi, xf);
		const Real w1 = lattice(xi+1, xf-1);
		// lerp values
		const auto fade = [](Real t) { return t*t*t*(t*(t*6 - 15) + 10); };
		return lerp(w0, w1, fade(xf));
	}

	template<typename Real, typename G> inline
	Real noise2(Real x, Real y, G lattice) {
		// setup
		const Real xi = floor(x);
		const Real yi = floor(y);
		const Real xf = x - xi;
		const Real yf = y - yi;
		// value = sample lattice
		const Real w00 = lattice(xi,   yi,   xf,   yf);
		const Real w10 = lattice(xi+1, yi,   xf-1, yf);
		const Real w01 = lattice(xi,   yi+1, xf,   yf-1);
		const Real w11 = lattice(xi+1, yi+1, xf-1, yf-1);
		// lerp values
		const auto fade = [](Real t) { return t*t*t*(t*(t*6 - 15) + 10); };
		return bilerp(w00, w10, w01, w11, fade(xf), fade(yf));
	}

	template<typename Real, typename G> inline 
	Real noise3(Real x, Real y, Real z, G lattice) {
		// setup
		const Real xi = floor(x);
		const Real yi = floor(y);
		const Real zi = floor(z);
		const Real xf = x - xi;
		const Real yf = y - yi;
		const Real zf = z - zi;
		// value = sample lattice
		const Real w000 = lattice(xi,   yi,   zi,   xf,   yf,   zf);
		const Real w100 = lattice(xi+1, yi,   zi,   xf-1, yf,   zf);
		const Real w010 = lattice(xi,   yi+1, zi,   xf,   yf-1, zf);
		const Real w110 = lattice(xi+1, yi+1, zi,   xf-1, yf-1, zf);
		const Real w001 = lattice(xi,   yi,   zi+1, xf,   yf,   zf-1);
		const Real w101 = lattice(xi+1, yi,   zi+1, xf-1, yf,   zf-1);
		const Real w011 = lattice(xi,   yi+1, zi+1, xf,   yf-1, zf-1);
		const Real w111 = lattice(xi+1, yi+1, zi+1, xf-1, yf-1, zf-1);
		// lerp values
		const auto fade = [](Real t) { return t*t*t*(t*(t*6 - 15) + 10); };
		return trilerp(
			w000, w100, w010, w110, 
			w001, w101, w011, w111, 
			fade(xf), fade(yf), fade(zf));
	}
	
	template<typename Real, typename G> inline 
	Real noise1(Real x, Real rep, G lattice) {
		// setup
		const Real xi = floor(x) % rep;
		const Real xf = fract(x);
		const Real xi1 = (xi + 1) % rep;
		// value = sample lattice
		const Real w0 = lattice(xi, xf);
		const Real w1 = lattice(xi1, xf-1);
		// lerp values
		const auto fade = [](Real t) { return t*t*t*(t*(t*6 - 15) + 10); };
		return lerp(w0, w1, fade(xf));
	}
	
	template<typename Real, typename G> inline
	Real noise2(Real x, Real y, Real repx, Real repy, G lattice) {
		// setup
		const Real xi = floor(x) % repx;
		const Real yi = floor(y) % repy;
		const Real xf = fract(x);
		const Real yf = fract(y);
		const Real xi1 = (xi+1) % repx;
		const Real yi1 = (yi+1) % repy;
		// value = sample lattice
		const Real w00 = lattice(xi,  yi,  xf,   yf);
		const Real w10 = lattice(xi1, yi,  xf-1, yf);
		const Real w01 = lattice(xi,  yi1, xf,   yf-1);
		const Real w11 = lattice(xi1, yi1, xf-1, yf-1);
		// lerp values
		const auto fade = [](Real t) { return t*t*t*(t*(t*6 - 15) + 10); };
		return bilerp(w00, w10, w01, w11, fade(xf), fade(yf));
	}

	struct PerlinLattice {
		static constexpr size_t  _NOISE_PERM_SIZE = 256;
		static constexpr uint8_t _NOISE_PERM[_NOISE_PERM_SIZE * 2] = {
			151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96,
			53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142,
			// Rest of noise permutation table
			8, 99, 37, 240, 21, 10, 23,
			190,  6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
			88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168,  68, 175, 74, 165, 71, 134, 139, 48, 27, 166,
			77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244,
			102, 143, 54,  65, 25, 63, 161,  1, 216, 80, 73, 209, 76, 132, 187, 208,  89, 18, 169, 200, 196,
			135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186,  3, 64, 52, 217, 226, 250, 124, 123,
			5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42,
			223, 183, 170, 213, 119, 248, 152,  2, 44, 154, 163,  70, 221, 153, 101, 155, 167,  43, 172, 9,
			129, 22, 39, 253,  19, 98, 108, 110, 79, 113, 224, 232, 178, 185,  112, 104, 218, 246, 97, 228,
			251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241,  81, 51, 145, 235, 249, 14, 239, 107,
			49, 192, 214,  31, 181, 199, 106, 157, 184,  84, 204, 176, 115, 121, 50, 45, 127,  4, 150, 254,
			138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180,
			151, 160, 137, 91, 90, 15,
			131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23,
			190,  6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
			88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168,  68, 175, 74, 165, 71, 134, 139, 48, 27, 166,
			77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244,
			102, 143, 54,  65, 25, 63, 161,  1, 216, 80, 73, 209, 76, 132, 187, 208,  89, 18, 169, 200, 196,
			135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186,  3, 64, 52, 217, 226, 250, 124, 123,
			5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42,
			223, 183, 170, 213, 119, 248, 152,  2, 44, 154, 163,  70, 221, 153, 101, 155, 167,  43, 172, 9,
			129, 22, 39, 253,  19, 98, 108, 110, 79, 113, 224, 232, 178, 185,  112, 104, 218, 246, 97, 228,
			251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241,  81, 51, 145, 235, 249, 14, 239, 107,
			49, 192, 214,  31, 181, 199, 106, 157, 184,  84, 204, 176, 115, 121, 50, 45, 127,  4, 150, 254,
			138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
		};

		template<typename Real>
		Real operator()(Real xi, Real xf) const {
			size_t ix = static_cast<ptrdiff_t>(xi) & (_NOISE_PERM_SIZE - 1);
			uint8_t h = _NOISE_PERM[ix] & 3;
			return (h & 1 ? -xf : xf);
		}

		template<typename Real>
		Real operator()(Real xi, Real yi, Real xf, Real yf) const {
			size_t ix = static_cast<ptrdiff_t>(xi) & (_NOISE_PERM_SIZE - 1);
			size_t iy = static_cast<ptrdiff_t>(yi) & (_NOISE_PERM_SIZE - 1);
			uint8_t h = _NOISE_PERM[_NOISE_PERM[ix] + iy] & 3;
			return ((h & 1) ? -xf : xf) + ((h & 2) ? -yf : yf);
		}

		template<typename Real>
		Real operator()(Real xi, Real yi, Real zi, Real xf, Real yf, Real zf) const {
			size_t ix = static_cast<ptrdiff_t>(xi) & (_NOISE_PERM_SIZE - 1);
			size_t iy = static_cast<ptrdiff_t>(yi) & (_NOISE_PERM_SIZE - 1);
			size_t iz = static_cast<ptrdiff_t>(zi) & (_NOISE_PERM_SIZE - 1);
			uint8_t h = _NOISE_PERM[ _NOISE_PERM[ _NOISE_PERM[ix] + iy ] + iz ] & 15;
			Real u = h < 8 || h == 12 || h == 13 ? xf : yf;
			Real v = h < 4 || h == 12 || h == 13 ? yf : zf;
			return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
		}
	};

	template<typename Real> inline
	Real perlin(Real x) {
		return noise1(x, PerlinLattice());
	}
	
	template<typename Real> inline
	Real perlin(Real x, Real y) {
		return noise2(x, y, PerlinLattice());
	}
	
	template<typename Real> inline
	Real perlin(Real x, Real y, Real z) {
		return noise2(x, y, z, PerlinLattice());
	}




	template<typename Real, typename Fn1, typename Fn2> inline
	Real vnoise1(Real x, Real rep, Fn1 hash, Fn2 fade) {
		// repeat lattice i
		Real i0 = floor( x ) % rep;
		Real i1 = (i0 + 1) % rep;
		// lerp
		return lerp( hash(i0), hash(i1), fade(fract(x)) );
	}

	template<typename Array2, typename Fhash21, typename F> inline
	auto vnoise2(Array2 x, Array2 rep, Fhash21 hash, F fade) {
		// repeat lattice i
		Array2 i0 = floor( x ) % rep;
		Array2 i1 = (i0 + 1) % rep;
		// bilerp
		Array2 uv = fade( fract( x ) );
		return lerp( 
			lerp( hash(i0),
				  hash(Array2{i1[0],i0[1]}), uv[0]), 
			lerp( hash(Array2{i0[0],i1[1]}), 
				  hash(Array2{i1[0],i1[1]}), uv[0]), uv[1] );
	}

	template<typename Array3, typename Fhash31, typename F> inline
	auto vnoise3(Array3 x, Array3 rep, Fhash31 hash, F fade) {
		// repeat lattice i
		Array3 i0 = floor( x ) % rep;
		Array3 i1 = (i0 + 1) % rep;
		// trilerp
		Array3 uvw = fade( fract( x ) );
		return lerp( 
			lerp( lerp( hash(i0),
						hash(Array3{i1[0],i0[1],i0[2]}), uvw[0]),
			      lerp( hash(Array3{i0[0],i1[1],i0[2]}),
					  	hash(Array3{i1[0],i1[1],i0[2]}), uvw[0]), uvw[1] ),
		    lerp( lerp( hash(Array3{i0[0],i0[1],i1[2]}),
						hash(Array3{i1[0],i0[1],i1[2]}), uvw[0]),
			      lerp( hash(Array3{i0[0],i1[1],i1[2]}),
						hash(i1), uvw[0]), uvw[1] ), uvw[2] );
	}

	template<typename Real> inline
	Real vnoise1(Real X, Real rep) {
		return vnoise1(X, rep, 
			[](Real i) { return fract(fract(i * 0.3183099f + 0.71f) * 17) * 2 - 1; },
			[](Real t) { return t*t*t*(t*(t*6 - 15) + 10); });
	}

	template<typename Array2> inline
	auto vnoise2(Array2 X, Array2 rep) {
		return vnoise2(X, rep,
			[](Array2 i){ i = fract(i * 0.3183099f + Array2{0.71f,0.113f}) * 17; return fract(i[0]*i[1]*(i[0]+i[1])) * 2 - 1; },
			[](Array2 t){ return t*t*t*(t*(t*6 - 15) + 10); });
	}

	template<typename Array3> inline
	auto vnoise3(Array3 X, Array3 rep) {
		/*{ "hash_function":{
			"part1":"from integer to floating",
			"part2":"add offset(avoid nearst Zero)",
			"part3":"get fraction of i, and scale"
		} }*/
		return vnoise3(X, rep,
			[](Array3 i) { i = fract(i * 0.3183099f + Array3{0.71f,0.113f,0.419f}) * 17; return fract(i[0]*i[1]*i[2]*(i[0]+i[1]+i[2])) * 2 - 1; },
			[](Array3 t) { return t*t*t*(t*(t*6 - 15) + 10); });
	}

	// vnoiseN ...



	// classic noise 2-demension, period is 289 times
	template<typename Vector2, typename Fn1, typename Fn2>
	auto cnoise2(Vector2 P, Fn1 permute, Fn2 fade) -> std::remove_cvref_t<decltype(P[0])> {
		using Real = std::remove_cvref_t<decltype(P[0])>;
		using Array4 = calculation::Array4<Real>;

		// period of gradient-vectors
		Array4 iP = floor( Array4{P[0],P[1],P[0],P[1]} ) + Array4{0,0,1,1};
			   iP = iP % 289;
		
		// gradient-vector = g(P,i)
		/*
		Vector2  Pi = floor(P);
		Vector2 nPi = Pi + 1;
				 Pi = Pi % 289;
				nPi = nPi % 289;
		Vector2 g00 = Vector2{ Pi[0], Pi[1]};
		Vector2 g10 = Vector2{nPi[0], Pi[1]};
		Vector2 g01 = Vector2{ Pi[0],nPi[1]};
		Vector2 g11 = Vector2{nPi[0],nPi[1]};
		
		Real i = permute(permute(g00[0]) + g00[1]);
		g00[0] = fract(i / 41) * 2 - 1;
		g00[1] = abs(g00[0]) - 0.5;
		g00[0] -= round(g00[0]);  
		g00 = normalize(g00);

		...
		*/
		Array4 ix = Array4{iP[0],iP[2],iP[0],iP[2]};
		Array4 iy = Array4{iP[1],iP[1],iP[3],iP[3]};
		Array4 i  = permute( permute( ix ) + iy );
		Array4 gx = fract( i / 41 ) * 2 - 1;
		Array4 gy = abs( gx ) - 0.5;// [-0.5, +0.5]
			   gx = gx - floor( gx + 0.5 );// [-0.5, +0.5]
		Vector2 g00 = Vector2{gx[0],gy[0]};
		Vector2 g10 = Vector2{gx[1],gy[1]};
		Vector2 g01 = Vector2{gx[2],gy[2]};
		Vector2 g11 = Vector2{gx[3],gy[3]};
		Array4 glen = sqrt( Array4{dot(g00,g00),dot(g10,g10),dot(g01,g01),dot(g11,g11)} );
		g00/=glen[0]; g10/=glen[1]; g01/=glen[2]; g11/=glen[3];
		
		// gradient-weight = dot(g(P,i), v(P,i))
		Vector2  Pf = fract(P);
		Vector2 nPf = Pf - 1;
		Real w00 = dot( g00, Vector2{ Pf[0], Pf[1]} );
		Real w10 = dot( g10, Vector2{nPf[0], Pf[1]} );
		Real w01 = dot( g01, Vector2{ Pf[0],nPf[1]} );
		Real w11 = dot( g11, Vector2{nPf[0],nPf[1]} );

		// bilerp gradient-weights
		Array4 uv = fade( Array4{Pf[0],Pf[1]} );
		Real result = lerp(
			lerp( w00, w10, uv[0] ), 
			lerp( w01, w11, uv[0] ), 
			uv[1] );
		return result * static_cast<Real>(2.1);
	}

	template<typename Vector2, typename Fn1, typename Fn2>
	auto cnoise2(Vector2 P, Vector2 rep, Fn1 permute, Fn2 fade) -> std::remove_cvref_t<decltype(P[0])> {
		using Real = std::remove_cvref_t<decltype(P[0])>;
		using Array4 = calculation::Array4<Real>;

		// period of gradient-vectors
		Array4 iP = floor( Array4{P[0],P[1],P[0],P[1]} ) + Array4{0,0,1,1};
			   iP = iP % Array4{rep[0],rep[1],rep[0],rep[1]};
			   iP = iP % 289;
		
		// gradient-vector = g(P,i)
		Array4 ix = Array4{iP[0],iP[2],iP[0],iP[2]};
		Array4 iy = Array4{iP[1],iP[1],iP[3],iP[3]};
		Array4 i  = permute( permute( ix ) + iy );
		Array4 gx = fract( i / 41 ) * 2 - 1;
		Array4 gy = abs( gx ) - 0.5;// [-0.5, +0.5]
			   gx = gx - floor( gx + 0.5 );// [-0.5, +0.5]
		Vector2 g00 = Vector2{gx[0],gy[0]};
		Vector2 g10 = Vector2{gx[1],gy[1]};
		Vector2 g01 = Vector2{gx[2],gy[2]};
		Vector2 g11 = Vector2{gx[3],gy[3]};
		Array4 glen = sqrt( Array4{dot(g00,g00),dot(g10,g10),dot(g01,g01),dot(g11,g11)} );
		g00/=glen[0]; g10/=glen[1]; g01/=glen[2]; g11/=glen[3];
		
		// gradient-weight = dot(g(P,i), v(P,i))
		Vector2  Pf = fract(P);
		Vector2 nPf = Pf - 1;
		Real w00 = dot( g00, Vector2{ Pf[0], Pf[1]} );
		Real w10 = dot( g10, Vector2{nPf[0], Pf[1]} );
		Real w01 = dot( g01, Vector2{ Pf[0],nPf[1]} );
		Real w11 = dot( g11, Vector2{nPf[0],nPf[1]} );

		// bilerp gradient-weights
		Array4 uv = fade( Array4{Pf[0],Pf[1]} );
		Real result = lerp(
			lerp( w00, w10, uv[0] ), 
			lerp( w01, w11, uv[0] ), 
			uv[1] );
		return result * static_cast<Real>(2.1);
	}

	// classic noise 2-demension, period is 289 times
	template<typename Vector2> inline
	auto cnoise2(Vector2 P) -> std::remove_cvref_t<decltype(P[0])> {
		return cnoise2(P, 
			[]<typename Ty>(Ty v){ return (v*v*34 + v) % 289; },
			[]<typename Ty>(Ty t){ return t*t*t*(t*(t*6 - 15) + 10); });
	}

	template<typename Vector2> inline
	auto cnoise2(Vector2 P, Vector2 rep) -> std::remove_cvref_t<decltype(P[0])> {
		return cnoise2(P, rep,
			[]<typename Ty>(Ty v){ return (v*v*34 + v) % 289; },
			[]<typename Ty>(Ty t){ return t*t*t*(t*(t*6 - 15) + 10); });
	}


	// classic noise 3-demension, rep is 289 example cnoise(P*289*Integer)
	template<typename Vector3, typename Fn1, typename Fn2>
	auto cnoise3(Vector3 P, Fn1 permute, Fn2 fade) -> std::remove_cvref_t<decltype(P[0])> {
		using Real = std::remove_cvref_t<decltype(P[0])>;
		using Array4 = calculation::Array4<Real>;

		// period of gradient-vectors
		Array4 iP0 = floor( Array4{P[0],P[1],P[2]} );
		Array4 iP1 = iP0 + 1;
			   iP0 = iP0 % 289;
			   iP1 = iP1 % 289;
		
		// gradient-vector = g(P,i)
		Array4 ix  = Array4{iP0[0],iP1[0],iP0[0],iP1[0]};
		Array4 iy  = Array4{iP0[1],iP0[1],iP1[1],iP1[1]};
		Real   iz0 = iP0[2];
		Real   iz1 = iP1[2];
		Array4 ixy  = permute( permute( ix ) + iy );
		Array4 ixy0 = permute( ixy + iz0 );
		Array4 ixy1 = permute( ixy + iz1 );
		Array4 gx0  = ixy0 / 7;
		Array4 gy0  = fract( floor( gx0 ) / 7 ) - 0.5;
			   gx0  = fract( gx0 );
		Array4 gz0  = - abs( gx0 ) - abs( gy0 ) + 0.5;
		Array4 sz0  = (gz0 <= 0);
			   gx0 -= sz0 * ( (gx0 >= 0) - 0.5 );
			   gy0 -= sz0 * ( (gy0 >= 0) - 0.5 );
		Array4 gx1  = ixy1 / 7;
		Array4 gy1  = fract( floor( gx1 ) / 7 ) - 0.5;
			   gx1  = fract( gx1 );
		Array4 gz1  = - abs( gx1 ) - abs( gy1 ) + 0.5;
		Array4 sz1  = (gz1 <= 0);
			   gx1 -= sz1 * ( (gx1 >= 0) - 0.5 );
			   gy1 -= sz1 * ( (gy1 >= 0) - 0.5 );
		Vector3 g000 = Vector3{gx0[0],gy0[0],gz0[0]};
		Vector3 g100 = Vector3{gx0[1],gy0[1],gz0[1]};
		Vector3 g010 = Vector3{gx0[2],gy0[2],gz0[2]};
		Vector3 g110 = Vector3{gx0[3],gy0[3],gz0[3]};
		Vector3 g001 = Vector3{gx1[0],gy1[0],gz1[0]};
		Vector3 g101 = Vector3{gx1[1],gy1[1],gz1[1]};
		Vector3 g011 = Vector3{gx1[2],gy1[2],gz1[2]};
		Vector3 g111 = Vector3{gx1[3],gy1[3],gz1[3]};
		Array4 glen0 = sqrt( Array4{dot(g000,g000),dot(g100,g100),dot(g010,g010),dot(g110,g110)} );
		Array4 glen1 = sqrt( Array4{dot(g001,g001),dot(g101,g101),dot(g011,g011),dot(g111,g111)} );
		g000/=glen0[0]; g100/=glen0[1]; g010/=glen0[2]; g110/=glen0[3];
		g001/=glen1[0]; g101/=glen1[1]; g011/=glen1[2]; g111/=glen1[3];

		// gradient-weight = dot(g(P,i), v(P,i))
		Vector3 v0 = fract( P );
		Vector3 v1 = v0 - 1;
		Real w000 = dot( g000, v0 );
		Real w100 = dot( g100, Vector3{v1[0],v0[1],v0[2]} );
		Real w010 = dot( g010, Vector3{v0[0],v1[1],v0[2]} );
		Real w110 = dot( g110, Vector3{v1[0],v1[1],v0[2]} );
		Real w001 = dot( g001, Vector3{v0[0],v0[1],v1[2]} );
		Real w101 = dot( g101, Vector3{v1[0],v0[1],v1[2]} );
		Real w011 = dot( g011, Vector3{v0[0],v1[1],v1[2]} );
		Real w111 = dot( g111, v1 );

		// trilerp gradient-weights
		Array4 uvw = fade( Array4{v0[0],v0[1],v0[2]} );
		Array4 result_xy = lerp( Array4{w000,w100,w010,w110}, Array4{w001,w101,w011,w111}, uvw[2] );
		Array4 result_x = lerp( Array4{result_xy[0],result_xy[1]}, Array4{result_xy[2],result_xy[3]}, uvw[1] );
		Real result = lerp( result_x[0], result_x[1], uvw[0] );
		return result * static_cast<Real>(2.0);
	}

	template<typename Vector3, typename Fn1, typename Fn2>
	auto cnoise3(Vector3 P, Vector3 rep, Fn1 permute, Fn2 fade) {
		using Real = std::remove_cvref_t<decltype(P[0])>;
		using Array4 = calculation::Array4<Real>;

		// period of gradient-vectors
		Array4 iP0 = floor( Array4{P[0],P[1],P[2]} ) % Array4{rep[0],rep[1],rep[2]};
		Array4 iP1 = (iP0 + 1) % Array4{rep[0],rep[1],rep[2]};
			   iP0 = iP0 % 289;
			   iP1 = iP1 % 289;
		
		// gradient-vector = g(P,i)
		Array4 ix  = Array4{iP0[0],iP1[0],iP0[0],iP1[0]};
		Array4 iy  = Array4{iP0[1],iP0[1],iP1[1],iP1[1]};
		Real   iz0 = iP0[2];
		Real   iz1 = iP1[2];
		Array4 ixy  = permute( permute( ix ) + iy );
		Array4 ixy0 = permute( ixy + iz0 );
		Array4 ixy1 = permute( ixy + iz1 );
		Array4 gx0  = ixy0 / 7;
		Array4 gy0  = fract( floor( gx0 ) / 7 ) - 0.5;
			   gx0  = fract( gx0 );
		Array4 gz0  = - abs( gx0 ) - abs( gy0 ) + 0.5;
		Array4 sz0  = (gz0 <= 0);
			   gx0 -= sz0 * ( (gx0 >= 0) - 0.5 );
			   gy0 -= sz0 * ( (gy0 >= 0) - 0.5 );
		Array4 gx1  = ixy1 / 7;
		Array4 gy1  = fract( floor( gx1 ) / 7 ) - 0.5;
			   gx1  = fract( gx1 );
		Array4 gz1  = - abs( gx1 ) - abs( gy1 ) + 0.5;
		Array4 sz1  = (gz1 <= 0);
			   gx1 -= sz1 * ( (gx1 >= 0) - 0.5 );
			   gy1 -= sz1 * ( (gy1 >= 0) - 0.5 );
		Vector3 g000 = Vector3{gx0[0],gy0[0],gz0[0]};
		Vector3 g100 = Vector3{gx0[1],gy0[1],gz0[1]};
		Vector3 g010 = Vector3{gx0[2],gy0[2],gz0[2]};
		Vector3 g110 = Vector3{gx0[3],gy0[3],gz0[3]};
		Vector3 g001 = Vector3{gx1[0],gy1[0],gz1[0]};
		Vector3 g101 = Vector3{gx1[1],gy1[1],gz1[1]};
		Vector3 g011 = Vector3{gx1[2],gy1[2],gz1[2]};
		Vector3 g111 = Vector3{gx1[3],gy1[3],gz1[3]};
		Array4 glen0 = sqrt( Array4{dot(g000,g000),dot(g100,g100),dot(g010,g010),dot(g110,g110)} );
		Array4 glen1 = sqrt( Array4{dot(g001,g001),dot(g101,g101),dot(g011,g011),dot(g111,g111)} );
		g000/=glen0[0]; g100/=glen0[1]; g010/=glen0[2]; g110/=glen0[3];
		g001/=glen1[0]; g101/=glen1[1]; g011/=glen1[2]; g111/=glen1[3];

		// gradient-weight = dot(g(P,i), v(P,i))
		Vector3 v0 = fract( P );
		Vector3 v1 = v0 - 1;
		Real w000 = dot( g000, v0 );
		Real w100 = dot( g100, Vector3{v1[0],v0[1],v0[2]} );
		Real w010 = dot( g010, Vector3{v0[0],v1[1],v0[2]} );
		Real w110 = dot( g110, Vector3{v1[0],v1[1],v0[2]} );
		Real w001 = dot( g001, Vector3{v0[0],v0[1],v1[2]} );
		Real w101 = dot( g101, Vector3{v1[0],v0[1],v1[2]} );
		Real w011 = dot( g011, Vector3{v0[0],v1[1],v1[2]} );
		Real w111 = dot( g111, v1 );

		// trilerp gradient-weights
		Array4 uvw = fade( Array4{v0[0],v0[1],v0[2]} );
		Array4 result_xy = lerp( Array4{w000,w100,w010,w110}, Array4{w001,w101,w011,w111}, uvw[2] );
		Array4 result_x = lerp( Array4{result_xy[0],result_xy[1]}, Array4{result_xy[2],result_xy[3]}, uvw[1] );
		Real result = lerp( result_x[0], result_x[1], uvw[0] );
		return result * static_cast<Real>(2.0);
	}

	// classic noise 3-demension, rep is 289 example cnoise(P*289*Integer)
	template<typename Vector3> inline
	auto cnoise3(Vector3 P) -> std::remove_cvref_t<decltype(P[0])> {
		return cnoise3(P, 
			[]<typename Ty>(Ty v){ return (v*v*34 + v) % 289; },
			[]<typename Ty>(Ty t){ return t*t*t*(t*(t*6 - 15) + 10); });
	}

	template<typename Vector3> inline
	auto cnoise3(Vector3 P, Vector3 rep) -> std::remove_cvref_t<decltype(P[0])> {
		return cnoise3(P, rep, 
			[]<typename Ty>(Ty v){ return (v*v*34 + v) % 289; },
			[]<typename Ty>(Ty t){ return t*t*t*(t*(t*6 - 15) + 10); });
	}
}// namespace fractional