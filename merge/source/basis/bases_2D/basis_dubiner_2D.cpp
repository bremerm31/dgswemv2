#include "../bases_2D.hpp"

namespace Basis {
	Array2D<double> Dubiner_2D::GetPhi(uint p, const std::vector<Point<2>>& coordinates) {
		uint n_pts = coordinates.size();

		Array2D<double> phi((p + 1)*(p + 2) / 2);

		std::vector<double> n1(n_pts);
		std::vector<double> n2(n_pts);

		for (uint i = 0; i < n_pts; i++) {
			n1[i] = 2 * (1 + coordinates[i][Z1]) / (1 - coordinates[i][Z2]) - 1;
			n2[i] = coordinates[i][Z2];

			if (coordinates[i][Z2] == 1) n1[i] = NAN; //singular point (-1,1);
		}

		uint m = 0;
		for (uint i = 0; i <= p; i++) {
			for (uint j = 0; j <= p - i; j++) {
				phi[m] = ComputePhi(i, j, n1, n2);
				m++;
			}
		}

		return phi;
	}

	Array3D<double> Dubiner_2D::GetDPhi(uint p, const std::vector<Point<2>>& coordinates) {
		uint n_pts = coordinates.size();

		Array3D<double> dphi((p + 1)*(p + 2) / 2);

		std::vector<double> n1(n_pts);
		std::vector<double> n2(n_pts);

		for (uint i = 0; i < n_pts; i++) {
			n1[i] = 2 * (1 + coordinates[i][Z1]) / (1 - coordinates[i][Z2]) - 1;
			n2[i] = coordinates[i][Z2];
		}

		uint m = 0;
		for (uint i = 0; i <= p; i++) {
			for (uint j = 0; j <= p - i; j++) {
				dphi[m] = ComputeDPhi(i, j, n1, n2);
				m++;
			}
		}

		return dphi;
	}

	std::pair<bool, Array2D<double>> Dubiner_2D::GetMinv(uint p) {
		std::pair<bool, Array2D<double>> m_inv(true, Array2D<double>(1)); //diagonal

		m_inv.second[0].reserve((p + 1)*(p + 2) / 2);
		for (uint i = 0; i <= p; i++) {
			for (uint j = 0; j <= p - i; j++) {
				m_inv.second[0].push_back((2 * i + 1)*(i + j + 1) / 2.0);
			}
		}

		return m_inv;
	}

	std::vector<double> Dubiner_2D::ComputePhi(uint p, uint q, const std::vector<double>& n1, const std::vector<double>& n2) {
		uint n_pts = n1.size(); //CHECK IF n1.size() = n2.size()

		std::vector<double> phi(n_pts);

		std::vector<double> psi_p(jacobi_polynomial(p, 0, 0, n1));
		std::vector<double> psi_pq(jacobi_polynomial(q, 2 * p + 1, 0, n2));

		for (uint i = 0; i < n_pts; i++) {
			phi[i] = psi_p[i] * pow((1 - n2[i]) / 2, p)*psi_pq[i];

			if (n1[i] == NAN) { //value of Dubiner polynomial at singular point (-1,1)
				if (p == 0) {
					phi[i] = q + 1;
				}
				else {
					phi[i] = 0;
				}
			}
		}

		return phi;
	}

	Array2D<double> Dubiner_2D::ComputeDPhi(uint p, uint q, const std::vector<double>& n1, const std::vector<double>& n2) {
		uint n_pts = n1.size(); //CHECK IF n1.size() = n2.size()

		Array2D<double> dphi_d(2);
		dphi_d[Z1].reserve(n_pts);
		dphi_d[Z2].reserve(n_pts);

		std::vector<double> psi_p(jacobi_polynomial(p, 0, 0, n1));
		std::vector<double> psi_pq(jacobi_polynomial(q, 2 * p + 1, 0, n2));

		std::vector<double> dpsi_p_dn1(jacobi_polynomial_derivative(p, 0, 0, n1));
		std::vector<double> dpsi_pq_dn2(jacobi_polynomial_derivative(q, 2 * p + 1, 0, n2));

		for (uint i = 0; i < n_pts; i++) {
			dphi_d[Z1].push_back((2.0 / (1.0 - n2[i])) * dpsi_p_dn1[i] * pow((1.0 - n2[i]) / 2.0, p)*psi_pq[i]);

			dphi_d[Z2].push_back(((1 + n1[i]) / (1.0 - n2[i])) * dpsi_p_dn1[i] * pow((1.0 - n2[i]) / 2.0, p)*psi_pq[i] +
				psi_p[i] * (pow((1.0 - n2[i]) / 2.0, p)*dpsi_pq_dn2[i] - (p / 2.0)*pow((1.0 - n2[i]) / 2.0, p - 1)*psi_pq[i]));
		}

		return dphi_d;
	}
}