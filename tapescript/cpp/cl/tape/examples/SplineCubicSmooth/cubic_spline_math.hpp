/*
Copyright (C) 2003-2015 CompatibL

This file is part of TapeScript, an open source library and tape encoding
standard for adjoint algorithmic differentiation (AAD), available from

http://github.com/compatibl/tapescript (source)
http://tapescript.org (documentation)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef cl_tape_impl_cubic_spline_math_hpp
#define cl_tape_impl_cubic_spline_math_hpp

namespace cl
{
	// Smooth cubic spline algorithm is implemented according to the Pollock's article.
	// (available at https://www.researchgate.net/publication/241171192_Smoothing_with_Cubic_Splines).
	// Variable names exactly follow the article, pp. 16-25.
	// auxiliaries stands for p, h, r, f, q.
	// uvw stands for u, v, w.
	namespace tapescript
	{

		// Initialise spline helper variables.
		template<typename vec, typename T>
		inline void calculate_helper_variables(const vec& x0, std::vector<std::vector<T>>& aux)
		{
			int n = x0.size();
			for (int i = 0; i < n - 1; i++)
				if (x0[i] >= x0[i + 1])
					cl::throw_("Incorrect data: x values must be sorted and distinct.");
			aux.resize(5);
			// p.
			aux[0].resize(n - 2);
			// h.
			aux[1].resize(n - 1);
			// r.
			aux[2].resize(n - 1);
			// f.
			aux[3].resize(n - 2);
			// q.
			aux[4].resize(n - 2);

			aux[1][0] = x0[1] - x0[0];
			aux[2][0] = 3 / aux[1][0];
			for (int i = 1; i < n - 1; ++i)
			{
				aux[1][i] = x0[i + 1] - x0[i];
				aux[2][i] = 3 / aux[1][i];
				aux[0][i - 1] = 2 * (aux[1][i - 1] + aux[1][i]);
				aux[3][i - 1] = -(aux[2][i - 1] + aux[2][i]);
			}
		}

		template<typename Array, typename T>
		inline void calculate_q(const Array& y0, std::vector<std::vector<T>>& aux)
		{
			int n = y0.size();
			for (int i = 1; i < n - 1; ++i)
				aux[4][i - 1] = 3 * (y0[i + 1] - y0[i]) / aux[1][i] - 3 * (y0[i] - y0[i - 1]) / aux[1][i - 1];
		}

		// Prepare system of linear equations: calculate all non-zero matrix elements.
		template<typename T>
		inline void prepare_linear_system(T mu, std::vector<std::vector<T>>& uvw, std::vector<std::vector<T>>& aux)
		{
			int n = aux[0].size() + 2;
			uvw.resize(3);
			uvw[0].resize(n - 2);
			uvw[1].resize(n - 3);
			uvw[2].resize(n - 4);
			std::vector<T>& u = uvw[0];
			std::vector<T>& v = uvw[1];
			std::vector<T>& w = uvw[2];
			std::vector<T>& p = aux[0];
			std::vector<T>& h = aux[1];
			std::vector<T>& r = aux[2];
			std::vector<T>& f = aux[3];
			for (int i = 0; i < n - 4; i++)
			{
				u[i] = mu * (std::pow(r[i], 2) + std::pow(f[i], 2) + std::pow(r[i + 1], 2)) + p[i];
				v[i] = mu * r[i + 1] * (f[i] + f[i + 1]) + h[i + 1];
				w[i] = mu * r[i + 1] * r[i + 2];
			}
			u[n - 4] = mu * (std::pow(r[n - 4], 2) + std::pow(f[n - 4], 2) + std::pow(r[n - 3], 2)) + p[n - 4];
			u[n - 3] = mu * (std::pow(r[n - 3], 2) + std::pow(f[n - 3], 2) + std::pow(r[n - 2], 2)) + p[n - 3];
			v[n - 4] = mu * r[n - 3] * (f[n - 4] + f[n - 3]) + h[n - 3];
		}

		// Solve system for spline.
		template<typename T>
		inline void solve_spline_five_diagonal_system(std::vector<std::vector<T>>& uvw, std::vector<T>& q_)
		{
			std::vector<T>& u = uvw[0];
			std::vector<T>& v = uvw[1];
			std::vector<T>& w = uvw[2];
			int n = uvw[0].size() + 2;
			// Factorization procedure.
			w[0] /= u[0];
			v[0] /= u[0];
			u[1] -= u[0] * std::pow(v[0], 2);
			v[1] = (v[1] - u[0] * v[0] * w[0]) / u[1];
			w[1] /= u[1];
			for (int i = 2; i < n - 4; i++)
			{
				T& ui = u[i];
				T& vi = v[i];
				T& wi = w[i];
				ui -= u[i - 2] * std::pow(w[i - 2], 2) + u[i - 1] * std::pow(v[i - 1], 2);
				vi = (vi - u[i - 1] * v[i - 1] * w[i - 1]) / ui;
				wi /= ui;
			}
			T& un = u[n - 4];
			T& unn = u[n - 3];
			T& wn = w[n - 5];
			un -= u[n - 6] * std::pow(w[n - 6], 2) + u[n - 5] * std::pow(v[n - 5], 2);
			v[n - 4] = (v[n - 4] - u[n - 5] * v[n - 5] * wn) / un;
			unn -= u[n - 5] * std::pow(wn, 2) + un * std::pow(v[n - 4], 2);
			// Forward substitution.
			q_[1] -= v[0] * q_[0];
			for (int i = 2; i < n - 2; ++i)
				q_[i] -= v[i - 1] * q_[i - 1] + w[i - 2] * q_[i - 2];
			for (int i = 0; i < n - 2; ++i)
				q_[i] /= u[i];
			// Back substitution.
			q_[n - 4] -= v[n - 4] * q_[n - 3];
			for (int i = n - 5; i >= 0; --i)
				q_[i] -= (v[i] * q_[i + 1] + w[i] * q_[i + 2]);
		}

		template<typename T>
		// Calculate spline B coefficients (solution of linear system of equations).
		inline void find_b_coefficients(std::vector<std::vector<T>>& uvw, std::vector<std::vector<T>>& aux)
		{
			std::vector<T>& q_ = aux[4];
			solve_spline_five_diagonal_system(uvw, q_);
		}

		// Calculate spline D coefficients.
		template<typename Array, typename T>
		inline std::vector<T> find_d_coefficients(T mu, const Array& y0, std::vector<std::vector<T>>& aux)
		{
			int n = y0.size();
			std::vector<T> result(n - 1);
			result[0] = y0[0] - mu * aux[2][0] * aux[4][0];
			result[1] = y0[1] - mu * (aux[3][0] * aux[4][0] + aux[2][1] * aux[4][1]);
			for (int i = 2; i < n - 2; ++i)
				result[i] = y0[i] - mu * (aux[2][i - 1] * aux[4][i - 2] + aux[3][i - 1] * aux[4][i - 1] + aux[2][i] * aux[4][i]);
			result[n - 2] = y0[n - 2] - mu * (aux[2][n - 3] * aux[4][n - 4] + aux[3][n - 3] * aux[4][n - 3]);
			return result;
		}

		// Calculate all spline coefficients.
		template<typename T>
		inline void find_all_coefficients(T mu, const std::vector<T>& D, const std::vector<std::vector<T>>& aux, std::vector<std::vector<T>>& coeffs)
		{
			// Compact variable name for number of data points.
			int n = D.size() + 1;
			coeffs.resize(3);
			for (int i = 0; i < 3; i++)
				coeffs[i].resize(n - 1);
			coeffs[0][0] = aux[4][0] / (3 * aux[1][0]);
			coeffs[1][0] = 0.0;
			coeffs[2][0] = (D[1] - D[0]) / aux[1][0] - aux[4][0] * aux[1][0] / 3;
			coeffs[0][1] = (aux[4][1] - aux[4][0]) / (3 * aux[1][1]);
			coeffs[1][1] = aux[4][0];
			coeffs[2][1] = aux[4][0] * aux[1][0] + coeffs[2][0];
			for (int i = 2; i < n - 2; i++)
			{
				coeffs[0][i] = (aux[4][i] - aux[4][i - 1]) / (3 * aux[1][i]);
				coeffs[1][i] = aux[4][i - 1];
				coeffs[2][i] = (aux[4][i - 1] + aux[4][i - 2]) * aux[1][i - 1] + coeffs[2][i - 1];
			}
			// Since q[n - 1] = 0.
			coeffs[0][n - 2] = (0 - aux[4][n - 3]) / (3 * aux[1][n - 2]);
			coeffs[1][n - 2] = aux[4][n - 3];
			coeffs[2][n - 2] = (aux[4][n - 3] + aux[4][n - 4]) * aux[1][n - 3] + coeffs[2][n - 3];
		}

		template <typename T>
		inline void mult_spline_special_factor(int const& n, T& lambda) {
			lambda *= std::pow(n - 1., 0.5);
		}

		template <typename Array, typename T>
		inline void find_estimated_variance(Array const& y0, T &std) {
			T variance = 0, diff, delta;
			T t = y0[1] - y0[0];
			int k = y0.size() - 1;
			for (int i = 1; i < k; i++)
			{
				delta = y0[i + 1] - y0[i];
				diff = i * delta - t;
				t += delta;
				variance += (diff * diff) / ((i + 1.0) * i);
			}
			std = variance / (k - 1);
		}

		template <typename Array, typename T>
		inline void spline_lambda(Array const& y0, T &lambda) {
			find_estimated_variance(y0, lambda);
			mult_spline_special_factor(y0.size(), lambda);
			lambda = 1. / (1. + lambda);
		}

		// Class of matrix of several nonzero diagonals.
		template <typename Array, typename T>
		class diagonal_matrix {
		public:
			std::vector<Array> matrix;
			int pos = 0;
			int row = 0;
			int col = 0;
			int size = 0;
			diagonal_matrix() {}
			diagonal_matrix(Array &fill, int p, int r, int c) {
				pos = p;
				row = r;
				col = c;
				size = fill.size();
				matrix.resize(size);
				for (int i = 0; i < size; ++i) {
					int k = std::min(row, col);
					if (i < pos) {
						k -= std::max(0, pos - i - std::max(0, row - col));
					}
					else {
						k -= std::max(0, i - pos - std::max(0, col - row));
					}
					matrix[i].resize(k, fill[i]);
				}
			}
			diagonal_matrix(std::vector<Array>& vec, int p, int r, int c) {
				matrix = vec;
				pos = p;
				row = r;
				col = c;
				size = vec.size();
			}
			void write() {
				std::cout << "\n\n";
				std::vector<Array> matr(row, Array(col));
				for (int i = 0; i < pos; ++i) {
					for (int j = 0; j < matrix[i].size(); ++j) {
						matr[pos - i + j][j] = matrix[i][j];
					}
				}
				for (int i = pos; i < matrix.size(); ++i) {
					for (int j = 0; j < matrix[i].size(); ++j) {
						matr[j][i - pos + j] = matrix[i][j];
					}
				}
				for (int i = 0; i < row; ++i) {
					for (int j = 0; j < col; ++j) {
						std::cout << matr[i][j] << " ";
					}
					std::cout << "\n";
				}
				std::cout << "\n";
			}
			void multiply(diagonal_matrix const& vec1, diagonal_matrix const& vec2) {
				row = vec1.row;
				col = vec2.col;
				int dim = std::min(row, col);
				int beg = std::max(1 - row, -vec1.pos - vec2.pos);
				int end = std::min(col, 1 + vec1.size - vec1.pos + vec2.size - vec2.pos);
				size = end - beg;
				matrix.resize(size);
				pos = -beg;
				for (int i = 0; i < std::min(pos, size); ++i) {
					matrix[i].resize(dim - std::max(pos - i - std::max(row - col, 0), 0), 0);
				}
				for (int i = std::max(pos, 0); i < size; ++i) {
					matrix[i].resize(dim - std::max(i - pos - std::max(col - row, 0), 0), 0);
				}
				for (int i = 0; i < vec1.size; i++) {
					for (int j = 0; j < vec2.size; j++) {
						int ind = i - vec1.pos + j - vec2.pos + pos;
						if (ind >= 0 && ind < size) {
							//std::cout << i << j << ind << "\n";
							int s = std::max(std::min(j - vec2.pos, vec1.pos - i), 0);
							int s1 = std::max(0, std::max(0, vec2.pos - j) - std::max(0, i - vec1.pos));
							int s2 = std::max(0, -std::max(0, vec2.pos - j) + std::max(0, i - vec1.pos));
							int last_k = std::min(matrix[ind].size() - s, std::min(vec1.matrix[i].size() - s1, vec2.matrix[j].size() - s2));
							//std::cout << s << s1 << s2 << last_k << "\n";
							for (int k = 0; k < last_k; k++) {
								matrix[ind][k + s] += vec1.matrix[i][k + s1] * vec2.matrix[j][k + s2];
							}
						}
					}
				}
			}
			void mult_vec(Array const& vec, Array& product) {
				product.resize(row, 0.);
				for (int i = 0; i < std::min(pos, size); ++i) {
					for (int j = 0; j < matrix[i].size(); ++j) {
						product[pos - i + j] += matrix[i][j] * vec[j];
					}
				}
				for (int i = std::max(0, pos); i < size; ++i) {
					for (int j = 0; j < matrix[i].size(); ++j) {
						product[j] += matrix[i][j] * vec[j + i - pos];
					}
				}
			}
			void remult_vec(Array& product) {
				Array vec;
				vec = product;
				mult_vec(vec, product);
			}
			void vec_mult(Array const& vec, Array& product) {
				product.resize(col, 0);
				for (int i = 0; i < std::min(pos, size); ++i) {
					for (int j = 0; j < matrix[i].size(); ++j) {
						product[j] += vec[j + pos - i] * matrix[i][j];
					}
				}
				for (int i = std::max(0, pos); i < size; ++i) {
					for (int j = 0; j < matrix[i].size(); ++j) {
						product[i - pos + j] += vec[j] * matrix[i][j];
					}
				}
			}
			void vec_remult(Array& product) {
				Array vec;
				vec = product;
				vec_mult(vec, product);
			}
			void add_matrix(diagonal_matrix &x) {
				for (int i = 0; i < x.size; ++i) {
					int ind = pos + i - x.pos;
					for (int j = 0; j < matrix[ind].size(); ++j) {
						matrix[ind][j] += x.matrix[i][j];
					}
				}
			}
			void scalar(T mu) {
				for (int i = 0; i < size; ++i) {
					for (int j = 0; j < matrix[i].size(); ++j) {
						matrix[i][j] *= mu;
					}
				}
			}
			void copy(diagonal_matrix const& mat) {
				size = mat.size;
				pos = mat.pos;
				row = mat.row;
				col = mat.col;
				matrix.resize(size);
				for (int i = 0; i < size; ++i) {
					matrix[i].resize(mat.matrix[i].size());
					for (int j = 0; j < matrix[i].size(); ++j) {
						matrix[i][j] = mat.matrix[i][j];
					}
				}
			}
			void mult(diagonal_matrix const& mat) {
				diagonal_matrix<Array, T> cop;
				cop.copy(*this);
				multiply(cop, mat);
			}
			void comult(diagonal_matrix const& mat) {
				diagonal_matrix<Array, T> cop;
				cop.copy(*this);
				multiply(mat, cop);
			}
			void product(std::vector<diagonal_matrix> const& matrices) {
				if (matrices.size() == 1) {
					copy(matrices[0]);
				}
				if (matrices.size() > 1) {
					multiply(matrices[0], matrices[1]);
					for (int i = 2; i < matrices.size(); ++i) {
						mult(matrices[i]);
					}
				}
			}
		};

		//helpful functions.
		template<typename Array, typename T>
		inline void scalar_product(Array const& first_vec, Array const& second_vec, T &mu) {
			mu = 0.;
			for (int i = 0; i < first_vec.size(); ++i) {
				mu += first_vec[i] * second_vec[i];
			}
		}

		template<typename Array>
		inline Array spline_special_shift(Array const& q) {
			Array S(q.size() - 1);
			for (int i = 0; i < S.size(); ++i) {
				S[i] = q[i] - q[i + 1];
			}
			return S;
		}

		template<typename Array>
		inline Array spline_special_coshift(Array const& q) {
			int n = q.size();
			Array S(n + 1);
			S[0] = q[0];
			for (int i = 1; i < n; ++i) {
				S[i] = q[i] - q[i - 1];
			}
			S[n] = -q[n - 1];
			return S;
		}

		template<typename Array>
		inline Array spline_special_mat_vec_mult(Array const& Qdiag, Array const& q) {
			Array S;
			int n = q.size();
			S.resize(n + 2);
			S[0] = Qdiag[0] * q[0];
			S[1] = Qdiag[1] * (q[1] - q[0]) - Qdiag[0] * q[0];
			for (int i = 2; i < n; ++i) {
				S[i] = (Qdiag[i - 1] * (q[i - 2] - q[i - 1]) + Qdiag[i] * (q[i] - q[i - 1]));
			}
			S[n] = (Qdiag[n - 1] * (q[n - 2] - q[n - 1]) - Qdiag[n] * q[n - 1]);
			S[n + 1] = (Qdiag[n] * q[n - 1]);
			return S;
		}

		template<typename Array>
		inline Array spline_special_mat_vec_comult(Array const& Qdiag, Array const& q) {
			Array S;
			int n = q.size();
			S.resize(n - 2);
			for (int i = 0; i < n - 2; ++i) {
				S[i] = (Qdiag[i] * (q[i] - q[i + 1]) + Qdiag[i + 1] * (q[i + 2] - q[i + 1]));
			}
			return S;
		}

		//!!!!!!!!!!!!!!!!!!!!!!!!//
		template<typename Array, typename It>
		inline Array spline_special_matP_vec_mult(Array const& x0, It &b) {
			int n = x0.size();
			Array Pb(n - 2);
			Pb[0] = 2 * (x0[2] - x0[0]) * b[0] + (x0[2] - x0[1]) * b[1];
			for (int i = 1; i < n - 3; ++i) {
				Pb[i] = (x0[i + 1] - x0[i]) * b[i - 1] + 2 * (x0[i + 2] - x0[i]) * b[i] + (x0[i + 2] - x0[i + 1]) * b[i + 1];
			}
			Pb[n - 3] = 2 * (x0[n - 1] - x0[n - 3]) * b[n - 3] + (x0[n - 2] - x0[n - 3]) * b[n - 4];
			return Pb;
		}

		template<typename Array, typename T>
		inline std::vector<Array> spline_special_tmat_build(T mu, Array const &Qdiag, Array const& tQdiag, Array const& dtx) {
			std::vector<Array> y;
			y.resize(5);
			int n = dtx.size() + 1;
			y[2].resize(n - 2);
			for (int i = 0; i < n - 2; ++i) {
				y[2][i] = 2 * (mu * ((tQdiag[i] * (Qdiag[i] + Qdiag[i + 1]) + 2 * tQdiag[i + 1] * Qdiag[i + 1]) +
					(Qdiag[i] * (tQdiag[i] + tQdiag[i + 1]))) + dtx[i + 1] + dtx[i]);
			}
			y[1].resize(n - 3);
			for (int i = 0; i < n - 3; ++i) {
				y[1][i] = dtx[i + 1] - mu * (Qdiag[i + 1] * (2 * tQdiag[i + 1] + tQdiag[i] + tQdiag[i + 2]) + tQdiag[i + 1] * (2 * Qdiag[i + 1] + Qdiag[i] + Qdiag[i + 2]));
			}
			y[0].resize(n - 4);
			for (int i = 0; i < n - 4; ++i) {
				y[0][i] = mu * (tQdiag[i + 1] * Qdiag[i + 2] + Qdiag[i + 1] * tQdiag[i + 2]);
			}
			y[3] = y[1]; y[4] = y[0];
			return y;
		}

		template<typename Array>
		inline std::vector<Array> copy_mat_from_array(int const &n, Array const& uvwq) {
			typedef decltype(std::begin(uvwq)) It;
			std::vector<It> uvw = { std::begin(uvwq), std::begin(uvwq) + n - 2, std::begin(uvwq) + 2 * n - 5, std::begin(uvwq) + 3 * n - 9 };
			std::vector<Array> uvw_(3);
			for (int i = 0; i < 3; ++i) {
				uvw_[i].resize(n - 2 - i);
				for (int j = 0; j < uvw_[i].size(); ++j) {
					uvw_[i][j] = uvw[i][j];
				}
			}
			return uvw_;
		}

		template <typename Array, typename T>
		inline void add_spesial_product(T scalar, Array &target, Array const& x, Array const& y) {
			for (int i = 0; i < target.size(); ++i) {
				target[i] += scalar * x[i] * y[i];
			}
		}

		template <typename Array>
		inline void special_multiply_vec(Array const& vec, Array &target) {
			for (int i = 0; i < target.size(); ++i) {
				target[i] = vec[i] * target[i];
			}
		}

		// Solve system for spline.
		template<typename Array, typename T>
		inline void solve_spline_five_diagonal_system(std::vector<Array>& uvw, Array& q_)
		{
			Array& u = uvw[0];
			Array& v = uvw[1];
			Array& w = uvw[2];
			int n = uvw[0].size() + 2;
			// Factorization procedure.
			w[0] /= u[0];
			v[0] /= u[0];
			u[1] -= u[0] * std::pow(v[0], 2);
			v[1] = (v[1] - u[0] * v[0] * w[0]) / u[1];
			w[1] /= u[1];
			for (int i = 3; i <= n - 4; i++)
			{
				T& ui = u[i - 1];
				T& vi = v[i - 1];
				T& wi = w[i - 1];
				ui -= u[i - 3] * std::pow(w[i - 3], 2) + u[i - 2] * std::pow(v[i - 2], 2);
				vi = (vi - u[i - 2] * v[i - 2] * w[i - 2]) / ui;
				wi /= ui;
			}
			T& un = u[n - 4];
			T& unn = u[n - 3];
			T& wn = w[n - 5];
			un -= u[n - 6] * std::pow(w[n - 6], 2) + u[n - 5] * std::pow(v[n - 5], 2);
			v[n - 4] = (v[n - 4] - u[n - 5] * v[n - 5] * wn) / un;
			unn -= u[n - 5] * std::pow(wn, 2) + un * std::pow(v[n - 4], 2);
			// Forward substitution.
			q_[1] -= v[0] * q_[0];
			for (int i = 2; i < n - 2; ++i)
				q_[i] += (-v[i - 1] * q_[i - 1] - w[i - 2] * q_[i - 2]);
			for (int i = 0; i < n - 2; ++i)
				q_[i] /= u[i];
			// Back substitution.
			q_[n - 4] -= v[n - 4] * q_[n - 3];
			for (int i = n - 5; i >= 0; --i)
				q_[i] -= (v[i] * q_[i + 1] + w[i] * q_[i + 2]);
		}

		// Calculates b coeffs.
		template <typename Array, typename T>
		inline Array find_b_uvw(T lambda, Array const& x0, Array const& y0, std::vector<Array> &uvw) {
			int n = x0.size();
			double mu = 2 * (1 - lambda) / (3 * lambda);
			Array Qdiag(n - 1), h(n - 1), vec, q(n - 2), S;
			for (int i = 0; i < h.size(); ++i) {
				h[i] = x0[i + 1] - x0[i];
				Qdiag[i] = 3 / h[i];
			}
			for (int i = 0; i < n - 2; ++i) {
				q[i] = Qdiag[i] * (y0[i] - y0[i + 1]) + Qdiag[i + 1] * (y0[i + 2] - y0[i + 1]);
			}
			uvw.resize(3);
			uvw[0].resize(n - 2);
			for (int i = 0; i < n - 2; ++i) {
				uvw[0][i] = 2 * (mu * (Qdiag[i] * (Qdiag[i] + Qdiag[i + 1]) + Qdiag[i + 1] * Qdiag[i + 1]) + h[i + 1] + h[i]);
			}
			uvw[1].resize(n - 3);
			for (int i = 0; i < n - 3; ++i) {
				uvw[1][i] = h[i + 1] - mu * Qdiag[i + 1] * (2 * Qdiag[i + 1] + Qdiag[i] + Qdiag[i + 2]);
			}
			uvw[2].resize(n - 4);
			for (int i = 0; i < n - 4; ++i) {
				uvw[2][i] = mu * Qdiag[i + 1] * Qdiag[i + 2];
			}
			std::vector<Array> uvw_;
			uvw_.resize(3);
			for (int i = 0; i < 3; ++i) {
				uvw_[i].resize(uvw[i].size());
				for (int j = 0; j < uvw[i].size(); ++j) {
					uvw_[i][j] = uvw[i][j];
				}
			}
			solve_spline_five_diagonal_system<Array, T>(uvw_, q);
			return q;
		}

		// Find major helper data.
		template <typename Array, typename T>
		inline void find_small_data(T lambda, Array const& x0, Array const& y0, Array &uvwq) {
			int n = x0.size();
			T mu = 2 * (1 - lambda) / (3 * lambda);
			Array h(n - 1);
			for (int i = 0; i < h.size(); ++i) {
				h[i] = x0[i + 1] - x0[i];
			}
			uvwq.resize(5 * n - 12);
			typedef decltype(std::begin(uvwq)) It;
			std::vector<It> uvw = { std::begin(uvwq), std::begin(uvwq) + n - 2, std::begin(uvwq) + 2 * n - 5, std::begin(uvwq) + 3 * n - 9, std::begin(uvwq) + 4 * n - 10 };
			for (int i = 0; i < n - 1; i++) {
				uvw[3][i] = 3 / h[i];
			}
			auto Qdiag = uvw[3];
			for (int i = 0; i < n - 2; ++i) {
				uvw[0][i] = 2 * (mu * (Qdiag[i] * (Qdiag[i] + Qdiag[i + 1]) + Qdiag[i + 1] * Qdiag[i + 1]) + h[i + 1] + h[i]);
			}
			for (int i = 0; i < n - 3; ++i) {
				uvw[1][i] = h[i + 1] - mu * Qdiag[i + 1] * (2 * Qdiag[i + 1] + Qdiag[i] + Qdiag[i + 2]);
			}
			for (int i = 0; i < n - 4; ++i) {
				uvw[2][i] = mu * Qdiag[i + 1] * Qdiag[i + 2];
			}


			Array b(n - 2);
			for (int i = 0; i < n - 2; ++i) {
				b[i] = Qdiag[i] * (y0[i] - y0[i + 1]) + Qdiag[i + 1] * (y0[i + 2] - y0[i + 1]);
			}
			std::vector<Array> uvw_;
			uvw_.resize(3);
			for (int i = 0; i < 3; ++i) {
				uvw_[i].resize(n - 2 - i);
				for (int j = 0; j < uvw_[i].size(); ++j) {
					uvw_[i][j] = uvw[i][j];
				}
			}
			solve_spline_five_diagonal_system<Array, T>(uvw_, b);
			for (int i = 0; i < n - 2; ++i) {
				uvw[4][i] = b[i];
			}
		}

		// Find spline using helper variables.
		template <typename Iterator, typename Array, typename T>
		inline Array find_b_with_uvwq(T lambda, Array const& y0, std::vector<Iterator> const& uvw) {
			double mu = 2 * (1 - lambda) / (3 * lambda);
			int n = y0.size();
			Array b(n - 2);
			Iterator Qdiag = uvw[3];
			for (int i = 0; i < n - 2; ++i) {
				b[i] = Qdiag[i] * (y0[i] - y0[i + 1]) + Qdiag[i + 1] * (y0[i + 2] - y0[i + 1]);
			}
			std::vector<Array> uvw_;
			uvw_.resize(3);
			for (int i = 0; i < 3; ++i) {
				uvw_[i].resize(n - 2 - i);
				for (int j = 0; j < uvw_[i].size(); ++j) {
					uvw_[i][j] = uvw[i][j];
				}
			}
			solve_spline_five_diagonal_system<Array, T>(uvw_, b);
			return b;
		}

		// Find spline using helper variables.
		template <typename Array, typename T>
		inline Array find_spline_with_uvwq(T lambda, Array const& x0, Array const& y0, Array &uvwq) {
			double mu = 2 * (1 - lambda) / (3 * lambda);
			int n = x0.size();
			typedef decltype(std::begin(uvwq)) It; //It b = std::begin(uvwq) + 4 * n - 10;
			std::vector<It> uvw = { std::begin(uvwq), std::begin(uvwq) + n - 2, std::begin(uvwq) + 2 * n - 5, std::begin(uvwq) + 3 * n - 9 };
			Array b = find_b_with_uvwq(lambda, y0, uvw);
			Array S;
			It Qdiag = std::begin(uvwq) + 3 * n - 9;
			S.resize(n);
			S[0] = y0[0] - mu * Qdiag[0] * b[0];
			S[1] = y0[1] - mu * (Qdiag[1] * (b[1] - b[0]) - Qdiag[0] * b[0]);
			for (int i = 2; i < n - 2; ++i) {
				S[i] = y0[i] - mu * (Qdiag[i - 1] * (b[i - 2] - b[i - 1]) + Qdiag[i] * (b[i] - b[i - 1]));
			}
			S[n - 2] = y0[n - 2] - mu * (Qdiag[n - 3] * (b[n - 4] - b[n - 3]) - Qdiag[n - 2] * b[n - 3]);
			S[n - 1] = y0[n - 1] - mu * (Qdiag[n - 2] * b[n - 3]);
			return S;
		}

		// Calculates spline.
		template <typename Array, typename T>
		inline Array find_spline(T lambda, Array const& x0, Array const& y0) {
			int n = x0.size();
			double mu = 2 * (1 - lambda) / (3 * lambda);
			Array Qdiag(n - 1), h(n - 1), vec, q(n - 2), S;
			for (int i = 0; i < h.size(); ++i) {
				h[i] = x0[i + 1] - x0[i];
			}
			for (int i = 0; i < n - 1; i++) {
				Qdiag[i] = 3 / h[i];
			}
			for (int i = 0; i < n - 2; ++i) {
				q[i] = Qdiag[i] * (y0[i] - y0[i + 1]) + Qdiag[i + 1] * (y0[i + 2] - y0[i + 1]);
			}
			std::vector<Array> uvw(3);
			uvw[0].resize(n - 2);
			for (int i = 0; i < n - 2; ++i) {
				uvw[0][i] = 2 * (mu * (Qdiag[i] * (Qdiag[i] + Qdiag[i + 1]) + Qdiag[i + 1] * Qdiag[i + 1]) + h[i + 1] + h[i]);
			}
			uvw[1].resize(n - 3);
			for (int i = 0; i < n - 3; ++i) {
				uvw[1][i] = h[i + 1] - mu * Qdiag[i + 1] * (2 * Qdiag[i + 1] + Qdiag[i] + Qdiag[i + 2]);
			}
			uvw[2].resize(n - 4);
			for (int i = 0; i < n - 4; ++i) {
				uvw[2][i] = mu * Qdiag[i + 1] * Qdiag[i + 2];
			}
			solve_spline_five_diagonal_system<Array, T>(uvw, q);
			S.resize(n);
			S[0] = y0[0] - mu * Qdiag[0] * q[0];
			S[1] = y0[1] - mu * (Qdiag[1] * (q[1] - q[0]) - Qdiag[0] * q[0]);
			for (int i = 2; i < n - 2; ++i) {
				S[i] = y0[i] - mu * (Qdiag[i - 1] * (q[i - 2] - q[i - 1]) + Qdiag[i] * (q[i] - q[i - 1]));
			}
			S[n - 2] = y0[n - 2] - mu * (Qdiag[n - 3] * (q[n - 4] - q[n - 3]) - Qdiag[n - 2] * q[n - 3]);
			S[n - 1] = y0[n - 1] - mu * (Qdiag[n - 2] * q[n - 3]);
			return S;
		}

		template <typename Array, typename T>
		inline Array a_find_spline(T lambda, Array const& x0, Array const& y0) {
			Array saved_info;
			find_small_data<Array, T>(lambda, x0, saved_info);
			return find_spline_with_uvwq<Array, T>(lambda, x0, y0, saved_info);
		}

		// For reverse y.
		template <typename Array, typename T>
		inline Array find_cospline(T lambda, Array const& x0, Array const& ts) {
			int n = x0.size();
			double mu = 2 * (1 - lambda) / (3 * lambda);
			Array Qdiag, xdif, dxdif, vec, q, S;
			Qdiag.resize(n - 1); xdif.resize(n - 3); dxdif.resize(n - 2);
			for (int i = 0; i < n - 1; i++) {
				Qdiag[i] = 3 / (x0[i + 1] - x0[i]);
			}
			for (int i = 0; i < xdif.size(); ++i) {
				xdif[i] = x0[i + 2] - x0[i + 1];
			}
			for (int i = 0; i < dxdif.size(); ++i) {
				dxdif[i] = 2 * (x0[i + 2] - x0[i]);
			}
			diagonal_matrix<Array, T>
				Q(std::vector<Array>{ Qdiag }, 0, n - 1, n - 1),
				first(Array{ 1, -1 }, 0, n - 2, n - 1),
				last(Array{ -1, 1 }, 1, n - 1, n - 2),
				second(Array{ 1, -1 }, 0, n - 1, n),
				third(Array{ -1, 1 }, 1, n, n - 1),
				P(std::vector<Array>{ xdif, dxdif, xdif}, 1, n - 2, n - 2),
				A, mat;

			A.product(std::vector<diagonal_matrix<Array, T>> { third, Q, last });
			A.scalar(mu);
			q = spline_special_mat_vec_comult(Qdiag, ts);
			mat.product(std::vector<diagonal_matrix<Array, T>> { first, Q, second });
			A.comult(mat);
			A.add_matrix(P);
			std::vector<Array> uvw{ A.matrix[2], A.matrix[1], A.matrix[0] };
			solve_spline_five_diagonal_system<Array, T>(uvw, q);
			S = spline_special_mat_vec_mult(Qdiag, q);
			for (int i = 0; i < n; ++i) {
				S[i] = ts[i] - mu * S[i];
			}
			return S;
		}

		// For reverse y.
		template <typename Array, typename T>
		inline Array find_cospline_with_uvwq(T lambda, Array const& x0, Array const& ts, Array const& uvwq) {
			int n = x0.size();
			double mu = 2 * (1 - lambda) / (3 * lambda);
			Array Qdiag, vec, q, S;
			Qdiag.resize(n - 1);
			for (int i = 0; i < n - 1; i++) {
				Qdiag[i] = 3 / (x0[i + 1] - x0[i]);
			}
			q = spline_special_mat_vec_comult(Qdiag, ts);
			std::vector<Array> uvw = copy_mat_from_array(n, uvwq);
			solve_spline_five_diagonal_system<Array, T>(uvw, q);
			S = spline_special_mat_vec_mult(Qdiag, q);
			for (int i = 0; i < n; ++i) {
				S[i] = ts[i] - mu * S[i];
			}
			return S;
		}

		// Calculates x derivatives of spline.
		template <typename Array, typename T>
		inline Array a_find_derivetives_x(T mu, Array const& x0, Array const& y0, Array const& tx, Array const& qq, std::vector<Array> &uvwq) {
			int n = x0.size();
			Array Qdiag, tQdiag, dif, ddif, xdif, dxdif, vec, q, tS;
			Qdiag.resize(n - 1); tQdiag.resize(n - 1); dif.resize(n - 3); ddif.resize(n - 2); xdif.resize(n - 3); dxdif.resize(n - 2);
			for (int i = 0; i < n - 1; i++) {
				Qdiag[i] = 3 / (x0[i + 1] - x0[i]);
				tQdiag[i] = -3 / std::pow(x0[i + 1] - x0[i], 2) * (tx[i + 1] - tx[i]);
			}
			for (int i = 0; i < dif.size(); ++i) {
				dif[i] = tx[i + 2] - tx[i + 1];
			}
			for (int i = 0; i < ddif.size(); ++i) {
				ddif[i] = 2 * (tx[i + 2] - tx[i]);
			}
			diagonal_matrix<Array, T> Q(std::vector<Array>{ Qdiag }, 0, n - 1, n - 1),
				tQ(std::vector<Array>{ tQdiag }, 0, n - 1, n - 1),
				first(Array{ 1, -1 }, 0, n - 2, n - 1),
				last(Array{ -1, 1 }, 1, n - 1, n - 2),
				second(Array{ 1, -1 }, 0, n - 1, n),
				third(Array{ -1, 1 }, 1, n, n - 1),
				P(std::vector<Array>{ dif, ddif, dif}, 1, n - 2, n - 2),
				y, mat1, mat2;

			y.product(std::vector<diagonal_matrix<Array, T>> { first, tQ, second });
			mat2.product(std::vector<diagonal_matrix<Array, T>> { third, Q, last });
			mat1.product(std::vector<diagonal_matrix<Array, T>> { third, tQ, last });
			y.mult_vec(y0, q);
			y.mult(mat2);
			for (int i = 0; i <= y.pos; ++i) {
				int k = 2 * y.pos - i;
				for (int j = 0; j < y.matrix[i].size(); ++j) {
					y.matrix[i][j] = mu * (y.matrix[i][j] + y.matrix[k][j]);
					y.matrix[k][j] = y.matrix[i][j];
				}
			}
			y.add_matrix(P);
			y.mult_vec(qq, vec);
			for (int i = 0; i < vec.size(); ++i) {
				q[i] -= vec[i];
			}
			mat1.mult_vec(qq, vec);
			solve_spline_five_diagonal_system<Array, T>(uvw, q);
			mat2.mult_vec(q, tS);
			for (int i = 0; i < n; ++i) {
				tS[i] = -mu  * (vec[i] + tS[i]);
			}
			return tS;
		}

		// Calculates spline estimation derivatives with respect to x input values.
		template <typename Array, typename T>
		inline Array find_derivetives_x(T mu, Array const& x0, Array const& y0, Array const& tx, Array const& qq, std::vector<Array> &uvw) {
			int n = x0.size();
			Array Qdiag, tQdiag, dtx, vec, q, tS;
			Qdiag.resize(n - 1); tQdiag.resize(n - 1); dtx.resize(n - 1); q.resize(n - 2);
			for (int i = 0; i < n - 1; i++) {
				Qdiag[i] = 3 / (x0[i + 1] - x0[i]);
				tQdiag[i] = -3 / std::pow(x0[i + 1] - x0[i], 2) * (tx[i + 1] - tx[i]);
				dtx[i] = tx[i + 1] - tx[i];
			}
			for (int i = 0; i < n - 2; ++i) {
				q[i] = tQdiag[i] * (y0[i] - y0[i + 1]) + tQdiag[i + 1] * (y0[i + 2] - y0[i + 1]);
			}
			std::vector<Array> y = spline_special_tmat_build<Array, T>(mu, Qdiag, tQdiag, dtx);
			diagonal_matrix<Array, T> z(y, 2, n - 2, n - 2);
			z.mult_vec(qq, vec);
			for (int i = 0; i < vec.size(); ++i) {
				q[i] -= vec[i];
			}
			vec = spline_special_mat_vec_mult<Array>(tQdiag, qq);
			solve_spline_five_diagonal_system<Array, T>(uvw, q);
			tS = spline_special_mat_vec_mult<Array>(Qdiag, q);
			for (int i = 0; i < n; ++i) {
				tS[i] = -mu  * (vec[i] + tS[i]);
			}
			return tS;
		}

		template <typename Array, typename T>
		inline Array solve_x_spline_derivative_with_uvwq(T lambda, Array const& x, Array const& y, Array const& tx, Array &uvwq) {
			T mu = 2 * (1 - lambda) / (3 * lambda);
			int n = x.size();
			typedef decltype(std::begin(uvwq)) It;
			std::vector<It> uvw = { std::begin(uvwq), std::begin(uvwq) + n - 2, std::begin(uvwq) + 2 * n - 5, std::begin(uvwq) + 3 * n - 9 };
			Array b = find_b_with_uvwq(lambda, y, uvw);
			std::vector<Array> uvw_(3);
			for (int i = 0; i < 3; ++i) {
				uvw_[i].resize(n - 2 - i);
				for (int j = 0; j < uvw_[i].size(); ++j) {
					uvw_[i][j] = uvw[i][j];
				}
			}
			return find_derivetives_x(mu, x, y, tx, b, uvw_);
		}

		template <typename Array, typename T>
		inline Array solve_x_spline_derivative(T lambda, Array const& x, Array const& y, Array const& tx) {
			double mu = 2 * (1 - lambda) / (3 * lambda);
			std::vector<Array> uvw;
			Array b = find_b_uvw<Array, T>(lambda, x, y, uvw);
			return a_find_derivetives_x<Array, T>(mu, x, y, tx, b, uvw);
		}

		// Calculates spline estimation derivatives with respect to x input values.
		template <typename Array, typename T>
		inline Array a_find_rev_derivetives_x(T mu, Array const& x0, Array const& y0, Array const& tS, Array const& b, std::vector<Array> uvw) {
			int n = x0.size();
			Array Qdiag, tQdiag, vec, q, tb, tx, left, right;
			Qdiag.resize(n - 1); tQdiag.resize(n - 1);
			double h = 0;
			for (int i = 0; i < n - 1; i++) {
				h = x0[i + 1] - x0[i];
				Qdiag[i] = 3 / h;
				tQdiag[i] = 3 / std::pow(h, 2);
			}
			diagonal_matrix<Array, T> Q(std::vector<Array>{ Qdiag }, 0, n - 1, n - 1),
				tQ(std::vector<Array>{ tQdiag }, 0, n - 1, n - 1),
				first(Array{ 1, -1 }, 0, n - 2, n - 1),
				last(Array{ -1, 1 }, 1, n - 1, n - 2),
				second(Array{ 1, -1 }, 0, n - 1, n),
				third(Array{ -1, 1 }, 1, n, n - 1),
				v, mat2, A;

			mat2.product(std::vector<diagonal_matrix<Array, T>> { third, Q, last });
			mat2.vec_mult(tS, tb);
			solve_spline_five_diagonal_system<Array, T>(uvw, tb);

			tx.resize(n - 1, 0);

			left = spline_special_shift(tS);
			special_multiply_vec(tQdiag, left);
			right = spline_special_coshift(b);
			add_spesial_product(1., tx, left, right);

			// Nabla A first part.                        
			right = spline_special_coshift(b);
			left = spline_special_mat_vec_mult(Qdiag, tb);
			left = spline_special_shift(left);
			special_multiply_vec(tQdiag, left);
			add_spesial_product(-mu, tx, left, right);

			// Nabla A second part.
			left = spline_special_coshift(tb);
			special_multiply_vec(tQdiag, left);
			right = spline_special_mat_vec_mult(Qdiag, b);
			right = spline_special_shift(right);
			add_spesial_product(-mu, tx, left, right);

			// Nabla Q for b.
			right = spline_special_shift(y0);
			add_spesial_product(1., tx, left, right);

			// Nabla P.
			for (int i = 1; i < n - 2; ++i) {
				tx[i] += 2 * (b[i - 1] * tb[i - 1] + b[i] * tb[i]) + b[i] * tb[i - 1] + b[i - 1] * tb[i];
			}
			tx[0] += 2 * b[0] * tb[0];
			tx[n - 2] += 2 * b[n - 3] * tb[n - 3];

			second.scalar(-mu);
			second.vec_remult(tx);
			return tx;
		}

		template <typename Array, typename T>
		inline Array find_rev_derivetives_x(T lambda, Array const& x0, Array const& y0, Array const& tS, Array const& b, Array const& uvwq) {
			T mu = 2 * (1 - lambda) / (3 * lambda);
			int n = x0.size();
			Array Qdiag, tQdiag, vec, q, tb, tx, left, right;
			Qdiag.resize(n - 1); tQdiag.resize(n - 1);
			double h = 0;
			for (int i = 0; i < n - 1; i++) {
				h = x0[i + 1] - x0[i];
				Qdiag[i] = 3 / h;
				tQdiag[i] = 3 / std::pow(h, 2);
			}
			typedef decltype(std::begin(uvwq)) It;
			std::vector<It> uvw = { std::begin(uvwq), std::begin(uvwq) + n - 2, std::begin(uvwq) + 2 * n - 5, std::begin(uvwq) + 3 * n - 9 };
			tb = find_b_with_uvwq(lambda, tS, uvw);

			tx.resize(n - 1, 0);
			// Nabla Q for S.
			left = spline_special_shift(tS);
			special_multiply_vec(tQdiag, left);
			right = spline_special_coshift(b);
			add_spesial_product(1., tx, left, right);

			// Nabla A first part.
			right = spline_special_coshift(b);
			left = spline_special_mat_vec_mult(Qdiag, tb);
			left = spline_special_shift(left);
			special_multiply_vec(tQdiag, left);
			add_spesial_product(-mu, tx, left, right);

			// Nabla A second part.
			left = spline_special_coshift(tb);
			special_multiply_vec(tQdiag, left);
			right = spline_special_mat_vec_mult(Qdiag, b);
			right = spline_special_shift(right);
			add_spesial_product(-mu, tx, left, right);
			// Nabla Q for b.
			right = spline_special_shift(y0);
			add_spesial_product(1., tx, left, right);
			// Nabla P.
			for (int i = 1; i < n - 2; ++i) {
				tx[i] += 2 * (b[i - 1] * tb[i - 1] + b[i] * tb[i]) + b[i] * tb[i - 1] + b[i - 1] * tb[i];
			}
			tx[0] += 2 * b[0] * tb[0];
			tx[n - 2] += 2 * b[n - 3] * tb[n - 3];

			q.resize(tx.size() + 1);
			q[0] = -mu * tx[0];
			for (int i = 1; i < tx.size(); ++i) {
				q[i] = mu * (tx[i - 1] - tx[i]);
			}
			q[tx.size()] = mu * tx[tx.size() - 1];
			return q;
		}

		template <typename Array, typename T>
		inline Array a_solve_x_spline_derivative(T lambda, Array const& x, Array const& y, Array const& tx) {
			T mu = 2 * (1 - lambda) / (3 * lambda);
			Array b = find_b<Array, T>(lambda, x, y);
			return b_find_derivetives_x<Array, T>(mu, x, y, tx, b);
		}

		// Find dmu / dy.
		template <typename Array, typename T>
		inline void mu_derivative(T &dmu, Array const& y0, Array const& ty) {
			int n = y0.size();
			dmu = -(y0[n - 1] - y0[0]) * (ty[n - 1] - ty[0]) / (n - 1.);
			for (int i = 1; i < n; ++i) {
				dmu += (y0[i] - y0[i - 1]) * (ty[i] - ty[i - 1]);
			}
			dmu *= 4. / (3. * (n - 2)) * std::pow(n - 1., 0.5);
		}

		// Find dS / dmu.
		template <typename Array, typename T>
		inline Array spline_mu_derivative(Array const& x0, Array const& y0, Array const& uvwq) {
			int n = y0.size();
			Array l(n - 1), Qdiag(n - 1);
			typedef decltype(std::begin(uvwq)) It;
			std::vector<It> uvw = { std::begin(uvwq), std::begin(uvwq) + n - 2, std::begin(uvwq) + 2 * n - 5, std::begin(uvwq) + 3 * n - 9, std::begin(uvwq) + 4 * n - 10 };
			It b = uvw[4];
			int start = 3 * n - 9;
			Array Pb(n - 2);
			Pb[0] = 2 * (x0[2] - x0[0]) * b[0] + (x0[2] - x0[1]) * b[1];
			for (int i = 1; i < n - 3; ++i) {
				Pb[i] = (x0[i + 1] - x0[i]) * b[i - 1] + 2 * (x0[i + 2] - x0[i]) * b[i] + (x0[i + 2] - x0[i + 1]) * b[i + 1];
			}
			Pb[n - 3] = 2 * (x0[n - 1] - x0[n - 3]) * b[n - 3] + (x0[n - 2] - x0[n - 3]) * b[n - 4];
			std::vector<Array> uvw_ = copy_mat_from_array(n, uvwq);
			solve_spline_five_diagonal_system<Array, T>(uvw_, Pb);
			for (int i = start; i < 4 * n - 10; ++i) {
				Qdiag[i - start] = uvwq[i];
			}
			return -spline_special_mat_vec_mult(Qdiag, Pb);
		}

		template <typename Array, typename T>
		inline Array solve_l_spline_derivative_with_uvwq(T lambda, Array const& x0, Array const& y0, Array const& ty, Array const& uvwq) {
			Array lambda_part = spline_mu_derivative<Array, T>(x0, y0, uvwq);
			T dmu;
			mu_derivative(dmu, y0, ty);
			return dmu * lambda_part;
		}

		template <typename Array, typename T>
		inline Array find_lambda_reverse(T const& lambda, Array const& x0, Array const& y, Array const& z, Array const& uvwq) {
			int n = y.size();
			typedef decltype(std::begin(uvwq)) It;
			std::vector<It> uvw = { std::begin(uvwq), std::begin(uvwq) + n - 2, std::begin(uvwq) + 2 * n - 5, std::begin(uvwq) + 3 * n - 9, std::begin(uvwq) + 4 * n - 10 };
			It b = uvw[4];
			Array Pb(n - 2), Qdiag(n - 1), S, u(n - 1);
			Pb[0] = 2 * (x0[2] - x0[0]) * b[0] + (x0[2] - x0[1]) * b[1];
			for (int i = 1; i < n - 3; ++i) {
				Pb[i] = (x0[i + 1] - x0[i]) * b[i - 1] + 2 * (x0[i + 2] - x0[i]) * b[i] + (x0[i + 2] - x0[i + 1]) * b[i + 1];
			}
			int start = 3 * n - 9;
			for (int i = start; i < 4 * n - 10; ++i) {
				Qdiag[i - start] = uvwq[i];
			}
			std::vector<Array> uvw_ = copy_mat_from_array(n, uvwq);
			solve_spline_five_diagonal_system<Array, T>(uvw_, Pb);
			S = spline_special_mat_vec_mult(Qdiag, Pb);
			T mu = 0.;
			scalar_product(z, S, mu);
			mu *= -4. / (3. * (n - 2));
			mult_spline_special_factor(n, mu);
			for (int i = 0; i < u.size(); ++i) {
				u[i] = y[i + 1] - y[i];
			}
			for (int i = 1; i < u.size(); ++i) {
				S[i] = u[i - 1] - u[i];
			}
			S[0] = 1. / (n - 1) * (y[n - 1] - y[0]) - u[0];
			S[n - 1] = -1. / (n - 1) * (y[n - 1] - y[0]) + u[n - 2];
			for (int i = 0; i < S.size(); ++i) {
				S[i] *= mu;
			}
			return S;
		}

		template <typename Array, typename T>
		inline Array solve_x_rev_spline_derivative(T lambda, Array const& x, Array const& y, Array const& tS) {
			T mu = 2 * (1 - lambda) / (3 * lambda);
			std::vector<Array> uvw;
			Array b = find_b_uvw<Array, T>(lambda, x, y, uvw);
			return a_find_rev_derivetives_x<Array, T>(mu, x, y, tS, b, uvw);
		}

		template <typename Array, typename T>
		inline Array solve_x_rev_spline_derivative_with_uvwq(T lambda, Array const& x, Array const& y, Array const& tS, Array const& uvwq) {
			T mu = 2 * (1 - lambda) / (3 * lambda);
			typedef decltype(std::begin(uvwq)) It;
			int n = tS.size();
			std::vector<It> uvw = { std::begin(uvwq), std::begin(uvwq) + n - 2, std::begin(uvwq) + 2 * n - 5, std::begin(uvwq) + 3 * n - 9 };
			Array b = find_b_with_uvwq(lambda, y, uvw);
			return find_rev_derivetives_x<Array, T>(lambda, x, y, tS, b, uvwq);
		}

		template <typename Array, typename T>
		inline Array spline_first_forward_with_uvwq(T lambda, Array const& x0, Array const& y0, Array const& x01, Array const& y01, Array& uvwq) {
			auto& yderiv = find_spline_with_uvwq(lambda, x0, y01, uvwq);//spline_vec(x0, y01);
			auto& xderiv = solve_x_spline_derivative_with_uvwq(lambda, x0, y0, x01, uvwq);//x_spline_deriv_forward(x0, y0, x01);
			auto& lderiv = solve_l_spline_derivative_with_uvwq(lambda, x0, y0, y01, uvwq);
			return yderiv + xderiv + lderiv;
		}

		template <typename Array, typename T>
		inline Array spline_first_reverse_with_uvwq(T lambda, Array const& x0, Array const& y0, Array const& z0, Array const& uvwq) {
			auto& pxy = find_cospline_with_uvwq(lambda, x0, z0, uvwq);
			auto& pxl = find_lambda_reverse(lambda, x0, y0, z0, uvwq);
			return pxy + pxl;
		}

		// Second order
		// Find d2S / dmu2
		template <typename Array, typename T>
		inline Array spline_mumu_derivative(T const& mu, Array const& x0, Array const& y0, Array const& uvwq) {
			int n = y0.size();
			//T mu = 2. / 3. * (1 - lambda) / lambda;
			Array l(n - 2), Qdiag(n - 1);
			typedef decltype(std::begin(uvwq)) It;
			std::vector<It> uvw = { std::begin(uvwq), std::begin(uvwq) + n - 2, std::begin(uvwq) + 2 * n - 5, std::begin(uvwq) + 3 * n - 9, std::begin(uvwq) + 4 * n - 10 };
			It b = uvw[4];
			for (int i = 0; i < n - 1; ++i) {
				Qdiag[i] = uvw[3][i];
			}

			Array Pb = spline_special_matP_vec_mult(x0, b);
			std::vector<Array> uvw_ = copy_mat_from_array(n, uvwq);
			solve_spline_five_diagonal_system<Array, T>(uvw_, Pb);
			Pb = spline_special_mat_vec_mult(Qdiag, Pb);
			Pb = spline_special_mat_vec_comult(Qdiag, Pb);
			uvw_ = copy_mat_from_array(n, uvwq);
			solve_spline_five_diagonal_system<Array, T>(uvw_, Pb);
			return spline_special_mat_vec_mult(Qdiag, Pb);
		}


		// Find spline using helper variables.
		template <typename Iterator, typename Array, typename T>
		inline Array find_b_with_uvwq_m(T mu, Array const& y0, std::vector<Iterator> const& uvw) {
			//double mu = 2 * (1 - lambda) / (3 * lambda);
			int n = y0.size();
			Array b(n - 2);
			Iterator Qdiag = uvw[3];
			for (int i = 0; i < n - 2; ++i) {
				b[i] = Qdiag[i] * (y0[i] - y0[i + 1]) + Qdiag[i + 1] * (y0[i + 2] - y0[i + 1]);
			}
			std::vector<Array> uvw_;
			uvw_.resize(3);
			for (int i = 0; i < 3; ++i) {
				uvw_[i].resize(n - 2 - i);
				for (int j = 0; j < uvw_[i].size(); ++j) {
					uvw_[i][j] = uvw[i][j];
				}
			}
			solve_spline_five_diagonal_system<Array, T>(uvw_, b);
			return b;
		}

		// Find major helper data.
		template <typename Array, typename T>
		inline void find_small_data_m(T mu, Array const& x0, Array const& y0, Array &uvwq) {
			int n = x0.size();
			//T mu = 2 * (1 - lambda) / (3 * lambda);
			Array h(n - 1);
			for (int i = 0; i < h.size(); ++i) {
				h[i] = x0[i + 1] - x0[i];
			}
			uvwq.resize(5 * n - 12);
			typedef decltype(std::begin(uvwq)) It;
			std::vector<It> uvw = { std::begin(uvwq), std::begin(uvwq) + n - 2, std::begin(uvwq) + 2 * n - 5, std::begin(uvwq) + 3 * n - 9, std::begin(uvwq) + 4 * n - 10 };
			for (int i = 0; i < n - 1; i++) {
				uvw[3][i] = 3 / h[i];
			}
			auto Qdiag = uvw[3];
			for (int i = 0; i < n - 2; ++i) {
				uvw[0][i] = 2 * (mu * (Qdiag[i] * (Qdiag[i] + Qdiag[i + 1]) + Qdiag[i + 1] * Qdiag[i + 1]) + h[i + 1] + h[i]);
			}
			for (int i = 0; i < n - 3; ++i) {
				uvw[1][i] = h[i + 1] - mu * Qdiag[i + 1] * (2 * Qdiag[i + 1] + Qdiag[i] + Qdiag[i + 2]);
			}
			for (int i = 0; i < n - 4; ++i) {
				uvw[2][i] = mu * Qdiag[i + 1] * Qdiag[i + 2];
			}


			Array b(n - 2);
			for (int i = 0; i < n - 2; ++i) {
				b[i] = Qdiag[i] * (y0[i] - y0[i + 1]) + Qdiag[i + 1] * (y0[i + 2] - y0[i + 1]);
			}
			std::vector<Array> uvw_;
			uvw_.resize(3);
			for (int i = 0; i < 3; ++i) {
				uvw_[i].resize(n - 2 - i);
				for (int j = 0; j < uvw_[i].size(); ++j) {
					uvw_[i][j] = uvw[i][j];
				}
			}
			solve_spline_five_diagonal_system<Array, T>(uvw_, b);
			for (int i = 0; i < n - 2; ++i) {
				uvw[4][i] = b[i];
			}
		}


		// Find spline using helper variables.
		template <typename Array, typename T>
		inline Array find_spline_with_uvwq_m(T mu, Array const& x0, Array const& y0, Array &uvwq) {
			//double mu = 2 * (1 - lambda) / (3 * lambda);
			int n = x0.size();
			typedef decltype(std::begin(uvwq)) It; //It b = std::begin(uvwq) + 4 * n - 10;
			std::vector<It> uvw = { std::begin(uvwq), std::begin(uvwq) + n - 2, std::begin(uvwq) + 2 * n - 5, std::begin(uvwq) + 3 * n - 9 };
			Array b = find_b_with_uvwq_m(mu, y0, uvw);
			Array S;
			It Qdiag = std::begin(uvwq) + 3 * n - 9;
			S.resize(n);
			S[0] = y0[0] - mu * Qdiag[0] * b[0];
			S[1] = y0[1] - mu * (Qdiag[1] * (b[1] - b[0]) - Qdiag[0] * b[0]);
			for (int i = 2; i < n - 2; ++i) {
				S[i] = y0[i] - mu * (Qdiag[i - 1] * (b[i - 2] - b[i - 1]) + Qdiag[i] * (b[i] - b[i - 1]));
			}
			S[n - 2] = y0[n - 2] - mu * (Qdiag[n - 3] * (b[n - 4] - b[n - 3]) - Qdiag[n - 2] * b[n - 3]);
			S[n - 1] = y0[n - 1] - mu * (Qdiag[n - 2] * b[n - 3]);
			return S;
		}

		// Find dS / dmu.
		template <typename Array, typename T>
		inline Array spline_mu_derivative_m(Array const& x0, Array const& y0, Array const& uvwq) {
			int n = y0.size();
			Array l(n - 1), Qdiag(n - 1);
			typedef decltype(std::begin(uvwq)) It;
			std::vector<It> uvw = { std::begin(uvwq), std::begin(uvwq) + n - 2, std::begin(uvwq) + 2 * n - 5, std::begin(uvwq) + 3 * n - 9, std::begin(uvwq) + 4 * n - 10 };
			It b = uvw[4];
			int start = 3 * n - 9;
			Array Pb = spline_special_matP_vec_mult(x0, b);
			std::vector<Array> uvw_ = copy_mat_from_array(n, uvwq);
			solve_spline_five_diagonal_system<Array, T>(uvw_, Pb);
			for (int i = start; i < 4 * n - 10; ++i) {
				Qdiag[i - start] = uvwq[i];
			}
			return -spline_special_mat_vec_mult(Qdiag, Pb);
		}

		// Find d2mu / dy2
		template <typename Array, typename T>
		inline void mumu_derivative(T &dmu, Array const& ty) {
			int n = ty.size();
			dmu = -(ty[n - 1] * ty[n - 1] + ty[0] * ty[0] - ty[n - 1] * ty[0]) / (n - 1.);
			dmu += ty[n - 1] * ty[n - 1] + ty[0] * ty[0] - ty[n - 1] * ty[n - 2];
			for (int i = 1; i < n - 1; ++i) {
				dmu += ty[i] * (2. * ty[i] - ty[i - 1]);
			}
			dmu *= 4. / (3. * (n - 2)) * std::pow(n - 1., 0.5);
		}

		// Find d2mu / dy2
		template <typename Array, typename T>
		inline Array spline_second_forward_with_uvwq(T lambda, Array const& x0, Array const& y0, Array const& y01, Array& uvwq) {
			int n = ty.size(), mu = 2 * (1 - lambda) / (3 * lambda), dmu, dmumu;
			mu_derivative(dmu, y0, y01);
			mumu_derivative(dmumu, y01);
			return dmu * dmu * spline_mumu_derivative(mu, x0, y0, uvwq)
				+ dmumu * spline_mu_derivative_m(x0, y0, uvwq)
				+ 2. * dmu * spline_mu_derivative_m(x0, y01, uvwq);
		}

		template <typename Array, typename T>
		inline Array debug_find(Array const& x0, Array const& y0) {
			T lambda;
			spline_lambda(y0, lambda);
			int n = ty.size();
			T mu = 2 * (1 - lambda) / (3 * lambda), dmu, dmumu;
			Array saved_info;
			find_small_data_m(mu, data.x, data.y, saved_info);
			return find_spline_with_uvwq_m(mu, x0, y0, uvwq);
		}
	}
}

#endif // cl_tape_impl_cubic_spline_math_hpp