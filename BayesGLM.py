#!/usr/bin/env python
import numpy as np
import pandas as pd
import getopt
import sys, os
from cmdstan import CmdStan
from multiprocessing import Pool
from itertools import izip

if __name__ == '__main__':
    n_processes = 2
    clustering_attribute = 'Cell_type'
    outfiles_path = ''

    optlist, args = getopt.gnu_getopt(sys.argv[1:], "hi:o:c:p:s:", ["help", "input=","output=","clus_attr=","processes=",'stanfolder='])

    if optlist== [] and args == []:
        print 'pystancef -i [INPUTFILE] -o [OUTPUTFOLDER] -c [CLUSTER_ATTRIBUTE] -p [THREADS] -s [CMDSTANFOLDER]'
        sys.exit()
    for opt, a in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ('-i', '--input'):
            input_path = a
        elif opt in ("-o", "--output"):
            outfiles_path = a
        elif opt in ('-c','--clus_attr'):
            clustering_attribute = str(a)
        elif opt in ('-p', '--processes'):
            n_processes = int(a)
        elif opt in ('-s', '--stanfolder'):
            stan_folder = str(a)
        else:
            assert False, "%s option is not supported" % opt
    
    # Determine the format of the input and read it
    input_format = input_path.split('.')[-1]
    if input_format == 'cef':
        # Load the cef file
        from Cef_tools import CEF_obj
        cef = CEF_obj()
        try:
            print 'Loading CEF'
            cef.readCEF(input_path)
        except:
            print 'Error in loading CEF file'

    elif input_format == 'loom':
        # Load the loom file
        import loompy
        ds = loompy.connect(input_path)
    else:
        print 'The format specified is not accepted. Please provide a .loom or .cef file'
        sys.exit()

    #Load Stan controller
    stan = CmdStan(stan_folder)
    # Move to Stan Folder 
    current_dir = os.path.realpath(os.curdir)
    os.chdir(stan_folder)
    
    # Check the existence of the compiled model otherwise compile
    if os.path.isfile(os.path.join('models', 'bayessianGLM')):
        print 'Model is already compiled.'
    else:
        print 'Compiled model was not found. Defining and compiling model.'
        model_code = '''
        # Bayesian Negative Binomial regression for single-cell RNA-seq

        data {
            int<lower=0> N;                 # number of outcomes
            int<lower=0> K;                 # number of predictors
            matrix<lower=0>[N,K] x;         # predictor matrix 
            int y[N];                       # outcomes
        }

        parameters {
            vector<lower=1>[K] beta;  # coefficients
            real<lower=0.001> r;  # overdispersion
        }

        model {	
            vector[N] mu;
            vector[N] rv;

            # priors
            r ~ cauchy(0, 1);
            beta ~ pareto(1, 1.5);

            # vectorize the scale
            for (n in 1:N) {
                rv[n] <- square(r + 1) - 1;
            }

            # regression
            mu <- x * (beta - 1) + 0.001;
            y ~ neg_binomial(mu ./ rv, 1 / rv[1]);
        }
            '''
        print model_code

        stan.compile("bayessianGLM", model_code)

        print 'Model is now saved for future use at %s' % os.path.join(stan_folder, "models")

    if input_format == 'cef':
        print 'Formating the model input'
        for i,v in izip(cef.col_attr_names, cef.col_attr_values):
            if i == clustering_attribute:
                predictor_list = v
            if 'total' in i.lower():
                total_molecules = [float(j) for j in v ]
                
        for i,v in izip(cef.row_attr_names, cef.row_attr_values):
            if 'gene' in i.lower():
                gene_names = v

        total = sum(total_molecules)/len(total_molecules)
        total_molecules_norm = [j/total for j in total_molecules]
        predictors = ['Size']
        for i in predictor_list:
            if i not in predictors:
                predictors.append(i)

        predictor_matrix = []
        for i, c_p in enumerate( predictor_list ):
            predictor_matrix.append( [total_molecules_norm[i]] + [float(c_p==p) for p in predictors[1:]] )

        N = len( predictor_matrix )
        K = len( predictors )
    else:
        ## TODO Implement a loom version
        ## TODO add the possibility to include custom predictor matrix to the model
        raise NotImplementedError

    
    def one_gene_model(name_gene, gene_vector, predictor_matrix, N, K):
        my_data = {'N': N, 'K': K, 'x': np.array(predictor_matrix), 'y': np.array(gene_vector)}
        # Create cmdstan controller object
        stan = CmdStan(stan_folder)
        # Sample and save traces
        stan.fit(model_name="bayessianGLM", data=my_data, output_name=name_gene,
        method="sample", additional_options=[], debug=True, parse_output=False,
        folder=os.path.split(input_path)[-1].split('.')[0]) 

    # Code to execute multiple treads
    print 'Passing the inputs to mutliple threads'
    try:
        p = Pool(processes=n_processes)
        for name_gene, gene_vector in izip(gene_names, cef.matrix):
            p.apply_async(one_gene_model,(name_gene, gene_vector, predictor_matrix, N, K))
        p.close()
        p.join()
    except KeyboardInterrupt:
        p.terminate()
        print "You cancelled the program!"
        sys.exit(1)