def write_eff_mass(out_name, eff_mass, eff_mass_err, xdata, osc = False):
    output = open(out_name, "w")
    output.write("#nt m m_err\n")
    for i in range(len(eff_mass)):
        if osc:
            output.write(str(xdata[i]) + " " + str(eff_mass[i]) + " " + str(eff_mass_err[i]) + '\n')
        else:
            output.write(str(xdata[i] + 0.5) + " " + str(eff_mass[i]) + " " + str(eff_mass_err[i]) + '\n')
    output.close()


def write_header(f, nstates, nstates_osc):
        f.write("#nt ")

        for i in range(nstates):
            f.write("m_no_" + str(i) + " ")
            f.write("m_no_" + str(i) + "_err ")

        for i in range(nstates_osc):
            f.write("m_osc_" + str(i) + " ")
            f.write("m_osc_" + str(i) + "_err ")

        for i in range(nstates):
            f.write("A_no_" + str(i) + " ")
            f.write("A_no_" + str(i) + "_err ")

        for i in range(nstates_osc):
            f.write("A_osc_" + str(i) + " ")
            f.write("A_osc_" + str(i) + "_err ")

        f.write("chi^2/d.o.f.")

        f.write(" AICc")

        f.write(" nt_up_lim\n")



def write_sample(out_name, ranges, sample_data, nstates, nstates_osc):
    with open(out_name, "w") as f:
        write_header(f, nstates = nstates,
                nstates_osc = nstates_osc)
        for i in range(len(ranges)):
            for res in sample_data[i]:
                write_line(f, ranges[i], res[0], res[1], res[2], res[3], nstates, nstates_osc)



def write_fit_mass(out_name, ranges, fit_res, fit_err, chi_dof = None, aicc = None,
        nstates = 1, nstates_osc = 0):

    with open(out_name, "w") as f:
    #one state ratio fit
        if len(fit_res[0]) == 1:
            f.write("#nt m m_err chi^2/d.o.f. AICc nt_up_lim\n")
            for i in range(len(ranges)):
                f.write(str(ranges[i][0]) + " "
                        + str(abs(fit_res[i][0])) + " " + str(fit_err[i][0])
                        + " " + str(chi_dof[i]) + " " + str(aicc[i]) + " " + str(ranges[i][1])
                        + '\n')

        #two state ratio fit
        elif len(fit_res[0]) == 3:
            f.write("#nt m1 m1_err m2 m2_err A2 A2_err chi^2/d.o.f. AICc nt_up_lim\n")
            for i in range(len(ranges)):
                f.write(str(ranges[i][0]) + " "
                        + str(abs(fit_res[i][1])) + " " + str(fit_err[i][1]) + " "
                        + str(abs(fit_res[i][2])) + " " + str(fit_err[i][2]) + " "
                        + str(fit_res[i][0]) + " " + str(fit_err[i][0])
                        + " " + str(chi_dof[i]) + " " + str(aicc[i])  + " " + str(ranges[i][1])
                        + '\n')

        #general fit
        else:
            write_header(f, nstates, nstates_osc)
            for i in range(len(ranges)):
                write_line(f, ranges[i], fit_res[i], fit_err[i], chi_dof[i], aicc[i], nstates,
                        nstates_osc)


def write_line(f, fit_range, fit_res, fit_err, chi_dof, aicc,
        nstates, nstates_osc):
        f.write(str(fit_range[0]) + " ")

        for no in range(nstates):
            f.write(str(abs(fit_res[2*no + 1])) + " " + str(fit_err[2*no + 1])
                    + " ")

        for osc in range(nstates_osc):
            osc += nstates
            f.write(str(abs(fit_res[2*osc + 1])) + " " + str(fit_err[2*osc + 1]) + " ")

        for no in range(nstates):
            f.write(str(fit_res[2*no]) + " " + str(fit_err[2*no]) + " ")

        for osc in range(nstates_osc):
            osc += nstates
            f.write(str(fit_res[2*osc]) + " " + str(fit_err[2*osc]) + " ")

        f.write(str(chi_dof) + " ")

        f.write(str(aicc) + " ")

        f.write(str(fit_range[1]) + '\n')
            

def write_pcov(out_name, fit_range, pcovs, nstates, nstates_osc):
    with open(out_name, "w") as f:
        f.write("#nt ")
        for i in range(nstates):
            f.write("p_A_no_" + str(i) + " ")
            f.write("p_m_no_" + str(i) + " ")


        for i in range(nstates_osc):
            f.write("p_A_osc_" + str(i) + " ")
            f.write("p_m_osc_" + str(i) + " ")
        f.write("nt_up_lim\n")
        for ntmin, pcov in enumerate(pcovs):
            for i in pcov:
                print(fit_range[ntmin][0], end = " ", file = f)
                for j in i:
                    print(j, end=" ", file = f)
                print(fit_range[ntmin][1], file = f)


def write_corr(out_name, xdata, ydata, ydata_err):
    output = open(out_name, "w")
    output.write("#nt G G_err\n")
    for i in range(len(ydata)):
        output.write(str(xdata[i]) + " " + str(ydata[i]) + " " + str(ydata_err[i]) + '\n')
    output.close()