import numpy as np
import os
import subprocess
from six.moves import zip, xrange

class CmdStan(object):
    def __init__(self, stan_folder):
        '''
        CmdStan Obect

        Arg:
        stan_folder (str):	the folder where cmdstn is 
        '''
        self.stan_folder = stan_folder

    def compile(self, model_name, model_code, debug=True):
        """
        Complie the model and save it for future as models/<model_name>.stan
        Rerunning with the same model name will cause the old model to be overritten
        
        Args:
            model_name (string):		Name of the model, without the '.stan' extension
            model_code (string):		Code for the model
            debug (bool):				If true, print the output from the Stan compiler
        """
        # Move to Stan Folder
        current_dir = os.path.realpath(os.curdir)
        os.chdir(self.stan_folder)

        # Colud check here if the file already exists

        # Create the model file
        with open("models/" + model_name + ".stan", "w") as model_file:
            model_file.write(model_code)

        # Make compile the model in C++
        process = subprocess.Popen(["make","models/" + model_name],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = process.communicate()
        if debug:
            print("*** Verbose output (call with debug=False to suppress) ***")
            print(output[0])
            print(output[1])
        if process.returncode != 0:
            raise RuntimeError("Model could not be compiled, return code: " + str(process.returncode))

        # Switch back to cmstan Folder
        os.chdir(current_dir)


    def save_Rdump(self, filename, data):
        '''
        Save a dictionary to R dumop format
        Args:
            filename (str)
            data (dict)
        '''
        with open(filename, "w") as data_file:
            for key, val in data.items():
                data_file.write(key + " <- ")
                if type(val) == np.ndarray:
                    if len(val.shape) == 1:
                        data_file.write("c(" + ",".join([str(x) for x in val])+ ")")
                    else:
                        temp = val.flatten('F')
                        data_file.write(("structure(c(" + ",".join([str(x) for x in val.flatten('F')])+ "),.Dim=c(%d,%d))") % val.shape)
                else:  # scalar
                    data_file.write(str(val))
                data_file.write("\n")
    
    def fit(self, model_name, data, output_name, method="sample", additional_options=[], debug=True, parse_output=False, folder=''):
        """
        Fit the model using the supplied data and save posterior samples (optionally loads them in memory)
        
        Args:
            model_name (string):		The name of the model (must be already compiled)
            data (dict): 				Input data as numpy arrays or simple scalars
            output_name (str):			A string a file will be saved as samples/rfolder/out_put_name.csv
            method (string): 			"sample", "optimize" or "variational"
            additional_options (list):  list containing options ready to pass to cmd stan through Popen, 
                                        the order has to respect the hierarchical system described in the manual 
            debug (bool):				If true, print verbose error messages instead of raising exceptions
            parse_output (bool):		If True runs self.parse_results(output_name) and return its output
            folder (str):				The project folder within samples where the file is saved (default, save directly in samples)
        """
        current_dir = os.path.realpath(os.curdir)
        os.chdir(self.stan_folder)
        # Make sure the output folder exists otherwise create it 
        if not os.path.isdir(os.path.join("samples", folder)):
            print('Creating the folder <cmdstanfolder>/samples/%s' % folder)
            os.makedirs(os.path.join("samples", folder)) 

        # Save the data in R dump() format, only temporarly
        data_file_name = "data/" + model_name + '_' +  output_name + ".data.R"
        self.save_Rdump(data_file_name, data)
        
        # Run the comand to start the cmdstan process
        process = subprocess.Popen(["models/" + model_name, method] +
        additional_options + ["data", "file=" + data_file_name,
        "output", "file=" + os.path.join( "samples", folder, output_name + ".csv")],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Deal with the process output
        output = process.communicate()
        if debug:
            print ("*** Verbose output (call with debug=False to suppress) ***")
            print (output[0])
            print (output[1])

        # Deal with possible errors
        if process.returncode != 0:
            if debug:
                print("ERROR: No samples were generated!")
                return {}
            else:
                raise RuntimeError("No samples were generated, return code: " + str(process.returncode))
        
        #Delete the R dump() file
        os.remove(data_file_name)

        if parse_output:
            out_ = self.parse_results("samples/" + output_name + ".csv")
        else:
            out_ = None

        # Go back to initial directory
        os.chdir(current_dir)

        return out_

    @staticmethod
    def parse_results(file_path):
        # Parse out the parameters in a dictionary
        params = []
        samples = []
        with open(file_path, "r") as infile:
            for line in infile:
                if line[0:2] == "lp":
                    params = line.strip().split(",")
                    continue
                if line[0] == "#" or line[0] == "\n":
                    continue
                samples.append([float(x) for x in line.split(",")])
        samples = np.array(samples)
        result = {}
        for ix in xrange(len(params)):
            result[params[ix]] = samples[:, ix]

        return result


def parse_stan_out(fname):
    class StanResults(object):
        pass
    res = StanResults()
    with open(fname) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            elif line.startswith("lp_"):
                header = np.array(line.rstrip("\n").split(","))
            elif line[0] in "0123456789":
                average = np.fromstring(line, sep=",")
                break
    res.header = header
    posteriors = np.loadtxt(fname, comments='#', delimiter=",", skiprows=32)

    variables_start = [0]
    variables_stop = []
    variables_shape = []
    shape_name = 1
    variables_names = [header[0].split(".")[0]]
    for i, name in enumerate(header):
        if name.split(".")[0] not in variables_names:
            nm, *shape_name = header[i - 1].split(".")
            variables_names.append(name.split(".")[0])
            variables_start.append(i)
            variables_stop.append(i)
            variables_shape.append([int(j) for j in shape_name][::-1])

    variables_stop.append(i + 1)
    nm, rr, cc = header[i].split(".")
    shape_name = [rr, cc]
    variables_shape.append([int(j) for j in shape_name][::-1])

    for i, k in enumerate(variables_names):
        val = np.array(posteriors[:, variables_start[i]:variables_stop[i]])
        avg = np.array(average[variables_start[i]:variables_stop[i]])
        head = np.array(header[variables_start[i]:variables_stop[i]])
        if variables_shape[i] == [] or len(variables_shape[i]) == 1:
            setattr(res, k, val)
            setattr(res, k + "_mean", avg)
            setattr(res, k + "_head", head)
        else:
            setattr(res, k.upper(), np.reshape(val, variables_shape[i] + [-1], order="F"))
            setattr(res, k.upper() + "_flat", val)
            setattr(res, k.upper() + "_mean", np.reshape(avg, variables_shape[i], order="F"))
            setattr(res, k.upper() + "_head", np.reshape(head, variables_shape[i], order="F"))
 
    return res
