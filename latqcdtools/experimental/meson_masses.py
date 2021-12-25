#!/usr/bin/env python3

# Meson masses corresponding to hisq connected correlators. Notation for spatial direction. Tables according to:
#   https://arxiv.org/pdf/1411.3018.pdf
#   http://inspirehep.net/record/1476934/files/PoS%28LATTICE%202015%29199.pdf
#   https://arxiv.org/pdf/1010.1216.pdf

m_scalar_osc = {#(M1) 1 0++ scalar
        'll' : 980, #a0
        'sl' : 1425, #K0*
        'cl' : 2318, #D0*
        'cs' : 2317.7, #Ds0*+-
        'cc' : 3414.75, #chis0
        }

m_scalar_no = {#(M1) gamma_3*gamma_5 0-+ scalar
        'll' : 139.57018, #pion+-
        'sl' : 493.677, #K+-
        'ss' : 686, #virt eta_s (https://arxiv.org/pdf/1111.1710.pdf)
        'cl' : 1869.61, #D+-
        'cs' : 1968.30, #Ds+-
        }

m_pdsclr_no = {#(M2) gamma_5 0-+ pseudo scalar
        'll' : 139.57018, #pion+-
        'sl' : 493.677, #K+-
        'ss' : 686, #virt eta_s (https://arxiv.org/pdf/1111.1710.pdf)
        'cl' : 1869.61, #D+-
        'cs' : 1968.30, #Ds+-
        'cc' : 2983.6, #eta_c
        }

m_pdsclr_osc = {#(M2) gamma_3 0-+ pseudo scalar
        }

m_vector_no = { #(M6, M7) gamma_1, gamma_2 1-- vector
        'll' : 775.26, #rho
        'sl' : 891.66, #K*
        'ss' : 1019.461, #phi
        'cl' : 2006.96, #D*
        'cs' : 2112.1, #Ds*+-
        'cc' : 3096.916, #J/psi
        }

m_vector_osc = { #(M6, M7) gamma_2*gamma_4, gamma_1*gamma_4 1+- tensor
        'll' : 1229.5, #b1
        'cc' : 3525.38, #hc
        }

m_axvctr_no = { #(M3, M4) gamma_1*gamma_3, gamma_2*gamma_3, 1-- vector
        'll' : 775.26, #rho
        'sl' : 1414, #K*
        'ss' : 1019.461, #phi
        'cl' : 2006.96, #D*
        'cs' : 2112.1, #Ds*+-
        'cc' : 3096.916, #J/psi
        }
        
m_axvctr_osc = { #(M3, M4) gamma_1*gamma_5, gamma_2*gamma_5 1++ axial vector
        'll' : 1230, #a1
        'sl' : 1272, #K1
        'ss' : 1426.4, #f1(1285)
        #'ss' : 1281.9, #f1
        'cl' : 2421.4, #D1
        'cs' : 2459.5, #Ds1+-
        'cc' : 3510.66, #Chic1
        }

m_name_scalar_osc = {#(M1) 1 0++ scalar
        'll' : '$a_0$',
        'sl' : '$K_0^*$',
        'cl' : '$D_0^*$',
        'cs' : '$D_{s0}^*\\pm$',
        'cc' : '$\\chi_{s0}$',
        }

m_name_scalar_no = {#(M1) gamma_3*gamma_5 0-+ scalar
        'll' : '$\\pi^{\\pm}$',
        'sl' : '$K^{\\pm}$',
        'ss' : '$\\eta_{s\\bar{s}}$',
        'cl' : '$D^{\\pm}$',
        'cs' : '$D_s^{pm}$',
        'cc' : '$\\eta_c$'
        }

m_name_pdsclr_no = {#(M2) gamma_5 0-+ pseudo scalar
        'll' : '$\\pi^{\\pm}$',
        'sl' : '$K^{\\pm}$',
        'ss' : '$\\eta_{s\\bar{s}}$',
        'cl' : '$D^{\\pm}$',
        'cs' : '$D_s^{pm}$',
        'cc' : '$\\eta_c$'
        }

m_name_pdsclr_osc = {#(M2) gamma_3 0-+ pseudo scalar
        }

m_name_vector_no = { #(M6, M7) gamma_1, gamma_2 1-- vector
        'll' : '$\\rho$',
        'sl' : '$K^*$',
        'ss' : '$\\phi$',
        'cl' : '$D^*$',
        'cs' : '$D_s^{*\\pm}$',
        'cc' : '$J/\\psi$',
        }

m_name_vector_osc = { #(M6, M7) gamma_2*gamma_4, gamma_1*gamma_4 1+- tensor
        'll' : '$b_1$',
        'cc' : '$h_c$',
        }

m_name_axvctr_no = { #(M3, M4) gamma_1*gamma_3, gamma_2*gamma_3, 1-- vector
        'll' : '$\\rho$',
        'sl' : '$K^*$',
        'ss' : '$\\phi$',
        'cl' : '$D^*$',
        'cs' : '$D_s^{*\\pm}$',
        'cc' : '$J/\\psi$',
        }
        
m_name_axvctr_osc = { #(M3, M4) gamma_1*gamma_5, gamma_2*gamma_5 1++ axial vector
        'll' : '$a_1$',
        'sl' : '$K_1$',
        'ss' : '$f_1$',
        'cl' : '$D_1$',
        'cs' : '$D_{s1}^\\pm$',
        'cc' : '$\\chi_{c1}$',
        }
