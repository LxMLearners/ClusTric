def write_tab_file(df, file_name, features, tps):
    with open(file_name, 'wt') as out:
        out.write("Total Times:\t" + str(tps) + '\n')
        out.write("Total Samples:\t" + str(len(features)) + '\n')
        out.write("Total Genes:\t" + str(len(df)) + '\n')

        for i in range(0, tps):
            fs = list(map(lambda x: str(i)+x, features)
                      ) if tps > 1 else features
            sub_feat = df[["Patient_ID"]+fs]

            out.write("Time\t" + str(i) + '\n')
            out.write("ID\tNAME\t")
            for l in range(len(features)):
                if l != len(features)-1:
                    out.write("S-" + str(l) + '\t')
                else:
                    out.write("S-" + str(l) + '\n')

            for index, row in sub_feat.iterrows():

                t_string = "{}\t"*(len(fs)+2)
                tupl = (index, "G-" + str(index)) + \
                    tuple(map(lambda x: row[x], fs))
                line = t_string.format(*tupl)
                out.write(line + '\n')

    out.close()
