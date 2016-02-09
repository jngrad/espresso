 /*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/ 
//comentario siguiente
#include "statistics_correlation_tcl.hpp"
#include "statistics_correlation.hpp"
#include "statistics_observable_tcl.hpp"
#include "statistics_observable.hpp"
#include "particle_data.hpp"
#include "parser.hpp"
#include "integrate.hpp" // needed to access the sim_time global variable

//  printf("\nOK here\n\n"); fflush(stdout);


/* forward declarations */
int tclcommand_print_correlation_time(Correlation& corr, Tcl_Interp* interp);
int tclcommand_print_average_errorbars(Correlation& corr, Tcl_Interp* interp); 
int tclcommand_analyze_parse_correlation(Tcl_Interp* interp, int argc, char** argv);
int tclcommand_correlation_parse_autoupdate(Tcl_Interp* interp, int no, int argc, char** argv);
int tclcommand_correlation_parse_print(Tcl_Interp* interp, int no, int argc, char** argv);
int tclcommand_correlation_parse_corr(Tcl_Interp* interp, int no, int argc, char** argv);
int double_correlation_print_spherically_averaged_sf(Correlation& corr, Tcl_Interp* interp);

/* General purpose and error reporting functions
************************************************
*/

// FIXME this is just a non-informative escape
int tclcommand_correlation_print_usage(Tcl_Interp* interp) {
  Tcl_AppendResult(interp, "You don't know how to use correlation.", (char *)NULL);
  return TCL_ERROR;
}


// FIXME this function is not used anywhere
//int observable_usage(Tcl_Interp* interp) {
//  Tcl_AppendResult(interp, "Usage: analyze correlation [ particle_velocities | com_velocity | particle_positions ", (char *) NULL);
//#ifdef ELECTROSTATICS
//  Tcl_AppendResult(interp, "particle_currents | currents ] ", (char *) NULL);
//#endif
//  Tcl_AppendResult(interp, "] [ first_obs | second_obs ] " , (char *)NULL);
//  Tcl_AppendResult(interp, "   type [ int ... int | * ]  \n" , (char *)NULL);
//  Tcl_AppendResult(interp, "or id [ int ... int ]  \n" , (char *)NULL);
//  Tcl_AppendResult(interp, "or molecule int \n" , (char *)NULL);
//  return TCL_ERROR;
//}


int correlation_print_parameters(Correlation& corr, Tcl_Interp* interp) {
  char buffer[16 + TCL_INTEGER_SPACE ];
  sprintf(buffer, " %d } ", corr.autocorrelation);
  Tcl_AppendResult(interp, "{ autocorrelation ", buffer, (char *)NULL);
  sprintf(buffer, " %d } ", corr.finalized);
  Tcl_AppendResult(interp, "{ finalized ", buffer, (char *)NULL);
  sprintf(buffer, " %d } ", corr.hierarchy_depth);
  Tcl_AppendResult(interp, "{ hierarchy_depth ", buffer, (char *)NULL);
  sprintf(buffer, " %d } ", corr.tau_lin);
  Tcl_AppendResult(interp, "{ tau_lin ", buffer, (char *)NULL);
  sprintf(buffer, " %d } ", corr.dim_A);
  Tcl_AppendResult(interp, "{ dim_A ", buffer, (char *)NULL);
  sprintf(buffer, " %d } ", corr.dim_B);
  Tcl_AppendResult(interp, "{ dim_B ", buffer, (char *)NULL);
  sprintf(buffer, " %d } ", corr.dim_corr);
  Tcl_AppendResult(interp, "{ dim_corr ", buffer, (char *)NULL);
  sprintf(buffer, " %d } ", corr.t);
  Tcl_AppendResult(interp, "{ t ", buffer, (char *)NULL);
  Tcl_PrintDouble(interp, corr.dt, buffer);
  Tcl_AppendResult(interp, "{ dt ", buffer, (char *)NULL);
  Tcl_AppendResult(interp, " } ", (char *)NULL);
  Tcl_PrintDouble(interp, corr.tau_max, buffer);
  Tcl_AppendResult(interp, "{ tau_max ", buffer, (char *)NULL);
  Tcl_AppendResult(interp, " } ", (char *)NULL);
  sprintf(buffer, " %d } ", corr.update_frequency);
  Tcl_AppendResult(interp, "{ update_frequency ", buffer, (char *)NULL);
  sprintf(buffer, " %d } ", corr.window_distance);
  Tcl_AppendResult(interp, "{ compressA ", corr.compressA_name, (char *)NULL);
  Tcl_AppendResult(interp, " } ", (char *)NULL);
  Tcl_AppendResult(interp, "{ compressB ", corr.compressB_name, (char *)NULL);
  Tcl_AppendResult(interp, " } ", (char *)NULL);
  Tcl_AppendResult(interp, "{ cor_operation ", corr.corr_operation_name, (char *)NULL);
  Tcl_AppendResult(interp, " } ", (char *)NULL);
  return TCL_OK;
}

int correlation_print_parameters_all(Tcl_Interp* interp) {
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE + 4];
  Tcl_AppendResult(interp, " { ", (char *)NULL);
  for (auto iter=correlations.begin(); iter !=correlations.end(); iter++) {
    sprintf(buffer, " %d ", iter->first);
    Tcl_AppendResult(interp, "\n{ correlation ", buffer, (char *)NULL);
    correlation_print_parameters(*(iter->second), interp);
  }
  Tcl_AppendResult(interp, " \n} ", (char *)NULL);
  return TCL_OK;
}

int correlation_print_average1(Correlation& corr, Tcl_Interp* interp, int argc, char** argv) {
  char buffer[TCL_DOUBLE_SPACE];
  int err;
  if (corr.n_data < 1) {
    Tcl_AppendResult(interp, buffer, "Error in print average: No input data available", (char *)NULL);
    return TCL_ERROR;
  }
  if (argc == 0) {
    for (unsigned i=0; i< corr.dim_A; i++) {
      Tcl_PrintDouble(interp, corr.A_accumulated_average[i]/corr.n_data, buffer);
      Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
    }
  return TCL_OK;
  } else if (ARG0_IS_S("formatted")) {
    double* values = (double*) Utils::malloc(corr.dim_A*sizeof(double));
    for (unsigned i=0; i< corr.dim_A; i++) {
      values[i]=corr.A_accumulated_average[i]/corr.n_data;
    }
    int change=0;
    err=tclcommand_observable_print_formatted(interp, argc-1, argv+1, &change, corr.A_obs, values);
    free(values);
    return err;
  } else {
    Tcl_AppendResult(interp, buffer, "Error in print average: No input data available", (char *)NULL);
    return TCL_ERROR;
  }
}

int correlation_print_variance1(Correlation& corr, Tcl_Interp* interp) {
  char buffer[TCL_DOUBLE_SPACE];
  if (corr.n_data < 1) {
    Tcl_AppendResult(interp, buffer, "Error in print variance: No input data available", (char *)NULL);
    return TCL_ERROR;
  }
  for (unsigned i=0; i< corr.dim_A; i++) {
    Tcl_PrintDouble(interp, 
        corr.A_accumulated_variance[i]/corr.n_data
        - (corr.A_accumulated_average[i]/corr.n_data)*(corr.A_accumulated_average[i]/corr.n_data), buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
  return TCL_OK;
}

int tclcommand_print_correlation_time(Correlation& corr, Tcl_Interp* interp) {
  char buffer[TCL_DOUBLE_SPACE];
  double* correlation_time;
  if (corr.dim_A != corr.dim_corr) {
    Tcl_AppendResult(interp, buffer, "Error in print correlation_time: Only makes sense when the dimensions of  \
        the observables and correlation match (Isn't it?) ", (char *)NULL);
    return TCL_ERROR;
  }
  if (!corr.autocorrelation) {
    Tcl_AppendResult(interp, buffer, "Does this make sense for an non-autocorrelation function? ", (char *)NULL);
    return TCL_ERROR;
  }
  correlation_time = (double*) Utils::malloc(corr.dim_corr*sizeof(double));
  corr.get_correlation_time(correlation_time);
  
  for (unsigned j=0; j<corr.dim_corr; j++) {
     Tcl_PrintDouble(interp, correlation_time[j], buffer);
     Tcl_AppendResult(interp, buffer, " ",(char *)NULL);
  }
  
  free(correlation_time);
  return TCL_OK;
}

int tclcommand_print_average_errorbars(Correlation& corr, Tcl_Interp* interp) {
  char buffer[TCL_DOUBLE_SPACE];
  double* correlation_time;
  double variance;
  double errorbar;
  if (corr.dim_A != corr.dim_corr) {
    Tcl_AppendResult(interp, buffer, "Error in print average_errorbars: Only makes sense when the dimensions of  \
        the observables and correlation match (Isn't it?) ", (char *)NULL);
    return TCL_ERROR;
  }
  if (!corr.autocorrelation) {
    Tcl_AppendResult(interp, buffer, "Error in print_average_errorbars: Must be an autocorrelation ", (char *)NULL);
    return TCL_ERROR;
  }
  correlation_time = (double*) Utils::malloc(corr.dim_corr*sizeof(double));
  corr.get_correlation_time(correlation_time);
  
  for (unsigned j=0; j<corr.dim_corr; j++) {
    variance=(corr.A_accumulated_variance[j]/corr.n_data - corr.A_accumulated_average[j]*corr.A_accumulated_average[j]/corr.n_data/corr.n_data);
    errorbar=sqrt(variance*(correlation_time[j]/corr.dt / corr.n_data));
    Tcl_PrintDouble(interp, errorbar, buffer);
    Tcl_AppendResult(interp, buffer, " ",(char *)NULL);
  }
  
  free(correlation_time);
  return TCL_OK;
}


/* We can track several correlation functions at a time
*  identified by their ids.
*/
int tclcommand_correlation(ClientData data, Tcl_Interp* interp, int argc, char** argv) {
  int no;
  char buffer[TCL_INTEGER_SPACE];
  argc-=1;
  argv+=1;
  if (argc < 1)
    return correlation_print_parameters_all(interp);
  if (argc==1 && ARG0_IS_S("n_corr")) {
    sprintf(buffer,"%d ",correlations.size());
    Tcl_AppendResult(interp, buffer, (char *)NULL);
    return TCL_OK;
  }
  if (ARG0_IS_S("new")) {
    // Find free id
    unsigned i=0;
   while ( correlations.find(i)!=correlations.end() )
     i++;
    argc-=1;
    argv+=1;
    return tclcommand_correlation_parse_corr(interp, i, argc, argv);
  } else if (ARG0_IS_I(no)) {
    if ( correlations.find(no) == correlations.end() ) {
      sprintf(buffer,"%d \n",no);
      Tcl_AppendResult(interp, "Correlation with id ", buffer, (char *)NULL);
      Tcl_AppendResult(interp, "is not defined", (char *)NULL);
      return TCL_ERROR;
    }
    argc-=1;
    argv+=1;
    return tclcommand_correlation_parse_corr(interp, no, argc, argv);
  } else {
    return tclcommand_correlation_print_usage(interp);
  }
}

int tclcommand_correlation_parse_autoupdate(Tcl_Interp* interp, int no, int argc, char** argv) {
  if (argc > 0 ) {
    if (ARG0_IS_S("start")) {
//      if(correlations.at(no)->A_fun==&tcl_input || correlations.at(no)->is_from_file) {
//        Tcl_AppendResult(interp, "Can not use autoupdate for a correlation from a file or from tclinput\n", (char *)NULL);
//        return TCL_ERROR;
//      }
      if(correlations.at(no)->update_frequency > 0) {
        correlations_autoupdate = 1;
        correlations.at(no)->autoupdate=1;
        correlations.at(no)->last_update=sim_time;
        return TCL_OK;
      } else {
        Tcl_AppendResult(interp, "Could not start autoupdate: update frequency not set\n", (char *)NULL);
        return TCL_ERROR;
      }
    } else if (ARG0_IS_S("stop")) {
      correlations_autoupdate=0;
      correlations.at(no)->autoupdate=0;
      for (auto iter=correlations.begin(); iter!=correlations.end(); iter++) {
        if (iter->second->autoupdate)
          correlations_autoupdate=1;
      }
      return TCL_OK;
    }
  }
  Tcl_AppendResult(interp, "Usage: analyze correlation $n autoupdate [ start | stop ]\n", (char *)NULL);
  return TCL_ERROR;
}


int tclcommand_correlation_parse_print(Tcl_Interp* interp, int no, int argc, char** argv) {
  if(argc==0) {
      return double_correlation_print_correlation(*(correlations.at(no)), interp);
  } 
  if(argc<1) { 
    Tcl_AppendResult(interp, "usage: analyze <correlation_id> print [what]\n", (char *)NULL); 
    return TCL_ERROR; 
  }
  if (ARG0_IS_S("spherically_averaged_sf")) { 
    return double_correlation_print_spherically_averaged_sf(*(correlations.at(no)), interp);
  } else if (ARG0_IS_S("average1")) {
    return correlation_print_average1(*(correlations.at(no)), interp, argc-1, argv+1);
  } else if (ARG0_IS_S("variance1")) {
    return correlation_print_variance1(*(correlations.at(no)), interp);
  } else if (ARG0_IS_S("correlation_time")) {
    return tclcommand_print_correlation_time(*(correlations.at(no)), interp);
  } else if (ARG0_IS_S("average_errorbars")) {
    return tclcommand_print_average_errorbars(*(correlations.at(no)), interp);
  } 
    
  Tcl_AppendResult(interp, "Unknown option to analyze correlation print: ", argv[0], (char *)NULL);
  return TCL_ERROR;

}

int tclcommand_correlation_parse_corr(Tcl_Interp* interp, int no, int argc, char** argv) {
//  int(*compressA)  ( double* A1, double*A2, double* A_compressed, unsigned int dim_A ) = 0;
 // int(*compressB)  ( double* B1, double*B2, double* B_compressed, unsigned int dim_B ) = 0;
  void **args = (void**)Utils::malloc(sizeof(void*)); // arguments to be passed to the correlation
  char *compressA_name=NULL;
  char *compressB_name=NULL;
  char *corr_operation_name=NULL;
  int tau_lin = 1; 
  double delta_t = 0.0;
  double tau_max = 0;
//  int(*corr_operation)  ( double* A, unsigned int dim_A, double* B, unsigned int dim_B, double* C, unsigned int dim_corr ) = 0;
  unsigned int dim_corr = 0;
  int change; // how many tcl argmuents are "consumed" by the parsing of arguments
  int error;
  char *error_msg=NULL;
  int temp;
  int obs_A_id,obs_B_id,dim_A,dim_B;
//  tcl_input_data tcl_input_d;
  char buffer[TCL_INTEGER_SPACE+TCL_DOUBLE_SPACE+2];
//  int autocorrelation=1; // by default, we are doing autocorrelation


  // Check if ID is negative
  if ( no < 0 ) {
    Tcl_AppendResult(interp, "Correlation IDs must be positive", (char *)NULL);
    return TCL_ERROR;
  }
  // Does it exist
  if ( correlations.find(no) != correlations.end()) {
    if (argc > 0) {
      if (ARG0_IS_S("print")) {
          return tclcommand_correlation_parse_print(interp, no, argc-1, argv+1); 
      } else if (ARG0_IS_S("print_params")) {
          if (argc > 1) {
              Tcl_AppendResult(interp, "\ncorrelation $id print_params does not accept further parameters\n", (char *)NULL);
              return TCL_ERROR;
          }
          return correlation_print_parameters(*(correlations.at(no)), interp);
      } else if (ARG0_IS_S("write_to_file")) {
        if (argc <2) {
          Tcl_AppendResult(interp, "You must pass a filename as argument of write_to_file", (char *)NULL);
          return TCL_ERROR;
        }
        error = double_correlation_write_to_file(*(correlations.at(no)), argv[1]);
        if (error) {
          Tcl_AppendResult(interp, "error in write to file", (char *)NULL);
        } else {
          return TCL_OK;
        }
      }
      // Operations beyond this point cannot be performed for finalized correlations 
      if ( correlations.at(no)->finalized ) {
	sprintf(buffer,"%d",no);
       Tcl_AppendResult(interp, "Correlation ", buffer, (char *)NULL);
        Tcl_AppendResult(interp, " has already been finalized. Cannot perform ", argv[0], (char *)NULL);
        Tcl_AppendResult(interp, "\n", buffer, (char *)NULL);
    	return TCL_ERROR;
      }
      if (ARG0_IS_S("finalize")) {
        if(argc==1) {
          return correlations.at(no)->finalize();
        } else {
          Tcl_AppendResult(interp, "Usage: analyze <correlation_id> finalize\n", buffer, (char *)NULL);
	  return TCL_ERROR;
        }
//      } else if (ARG0_IS_S("update") && correlations.at(no)->A_fun==&tcl_input ) {
//        correlations.at(no)->A_args = &tcl_input_d;
//        tcl_input_d.interp=interp;
//        tcl_input_d.argv=argv+1;
//        tcl_input_d.argc=argc-1;
//        if (!correlations.at(no)->autocorrelation) {
//          correlations.at(no)->B_args = &tcl_input; 
//          tcl_input_d.interp=interp;
//          tcl_input_d.argv=argv+2;
//          tcl_input_d.argc=argc-2;
//        }
//        error = double_correlation_get_data(correlations.at(no));
//        if (error) {
//          Tcl_AppendResult(interp, "Error reading tclinput", (char *)NULL);
//        } else {
//          return TCL_OK;
//        }
      } else if (ARG0_IS_S("update")) {
	if (correlations.at(no)->autoupdate == 1) {
	  sprintf(buffer, "%d. ", no);
          Tcl_AppendResult(interp, "Cannot update correlation ", buffer, (char *)NULL);
          Tcl_AppendResult(interp, "It is already in the autoupdate mode.\n", (char *)NULL);
	  return TCL_ERROR;
	}
        error = correlations.at(no)->get_data();
        if (error) {
          Tcl_AppendResult(interp, double_correlation_get_data_errors[error], (char *)NULL);
	  return TCL_ERROR;
        } else {
          return TCL_OK;
        }
      } else if (ARG0_IS_S("autoupdate")) {
        return tclcommand_correlation_parse_autoupdate(interp, no, argc-1, argv+1);
//      } else if (ARG0_IS_S("update_from_file")) {
//        error = correlation_update_from_file(no);
//        if (error) {
//          Tcl_AppendResult(interp, "error in update_from_file", (char *)NULL);
//          return TCL_ERROR;
//        } else {
//          return TCL_OK;
//        }
      } else if (ARG0_IS_S("write_checkpoint_binary")) {
	if (correlations.at(no)->write_data_to_file( argv[1],true)){
          Tcl_AppendResult(interp, "Failed to write checkpoint.\n", (char *)NULL);
	  return TCL_ERROR;
	}
	return TCL_OK;
      } else if (ARG0_IS_S("write_checkpoint_ascii")) {
	if ( correlations.at(no)->write_data_to_file(argv[1],false)){
          Tcl_AppendResult(interp, "Failed to write checkpoint.\n", (char *)NULL);
	  return TCL_ERROR;
	}
	return TCL_OK;
      } else if (ARG0_IS_S("read_checkpoint_binary")) {
	if (correlations.at(no)->read_data_from_file(argv[1],true)){
          Tcl_AppendResult(interp, "Failed to read checkpoint.\n", (char *)NULL);
	  return TCL_ERROR;
	}
	return TCL_OK;
      } else if (ARG0_IS_S("read_checkpoint_ascii")) {
	if (correlations.at(no)->read_data_from_file(argv[1],false)){
          Tcl_AppendResult(interp, "Failed to read checkpoint.\n", (char *)NULL);
	  return TCL_ERROR;
	}
	return TCL_OK;
      } else {
        Tcl_AppendResult(interp, "Usage for an already existing correlation:", (char *)NULL);
        Tcl_AppendResult(interp, "analyze correlation $ID [ print | update ]", (char *)NULL);
        return TCL_ERROR;
      }
    } else {
      Tcl_AppendResult(interp, "Usage for an already existing correlation:", (char *)NULL);
      Tcl_AppendResult(interp, "analyze correlation $ID [ print | update ]", (char *)NULL);
      return TCL_ERROR;
    }

  } 
  // Enforce consecutive ids
  else if ( (no == 0) || correlations.find(no-1) != correlations.end()) {  
  
    //Tcl_AppendResult(interp, "Setting up a new correlation\n", (char *)NULL);
    // Else we must parse the other arguments and see if we can construct a fully
    // working correlation class instance from that.
    while (argc > 0) {
      if ( ARG0_IS_S("first_obs") || ARG0_IS_S_EXACT("obs1") ) {
        if (argc>1 && ARG1_IS_I(temp)) {
          if (observables.find(temp)==observables.end()) {
             Tcl_AppendResult(interp, "Error in correlation observable. The specified observable does not exist\n", (char *)NULL);
             return TCL_ERROR;
          }
          obs_A_id=temp;
	  dim_A=observables.at(obs_A_id)->n;
          change+=2; argv+=2; argc-=2;
        } else {
          tclcommand_correlation_print_usage(interp);
          return TCL_ERROR;
        }
      } else if ( ARG0_IS_S("second_obs") || ARG0_IS_S_EXACT("obs2") ) {
        if (argc>1 && ARG1_IS_I(temp)) {
          if (observables.find(temp)==observables.end()) {
             Tcl_AppendResult(interp, "Error in correlation observable. The specified observable does not exist\n", (char *)NULL);
             return TCL_ERROR;
          }
          obs_B_id=temp;
	  dim_B=observables.at(obs_B_id)->n;
          change+=2; argv+=2; argc-=2;
	    //autocorrelation=0;
        } else {
          tclcommand_correlation_print_usage(interp);
          return TCL_ERROR;
        }
      } else if ( ARG0_IS_S("corr_operation") ) {
        argc -= 1;
        argv += 1;
        if ( parse_corr_operation(interp, argc, argv, &change, &corr_operation_name, &dim_corr, dim_A, dim_B, args) ) 
          return TCL_ERROR; 
        argc -= change;
        argv += change;
      } else if ( ARG0_IS_S("tau_lin") ) {
          if ( argc < 2 || ! (ARG1_IS_I(tau_lin)) ) 
              Tcl_AppendResult(interp, "Usage: analyze correlation ... tau_lin $tau_lin", (char *)NULL); 
          else { 
              argc -= 2; 
              argv += 2; 
          }
	  // tau_lin is already set
      } else if ( ARG0_IS_S("tau_max") ) {
        if ( argc < 2 || !ARG1_IS_D(tau_max)) {
          Tcl_AppendResult(interp, "Usage: analyze correlation ... tau_max $tau_max\n", (char *)NULL);
	  return TCL_ERROR;
	} else { 
          argc -= 2;
          argv += 2;
        }
      } else if ( ARG0_IS_S("delta_t") || ARG0_IS_S("dt") ) {
        if ( delta_t > 0 ) {
          Tcl_AppendResult(interp, "Error setting correlation: attempt to set delta_t twice", (char *)NULL);
	  return TCL_ERROR;
        } else if ( argc < 2 || !ARG1_IS_D(delta_t)) {
          Tcl_AppendResult(interp, "Usage: analyze correlation ... delta_t $delta_t ", (char *)NULL);
	  return TCL_ERROR;
        } else { 
	  argc -= 2; 
	  argv += 2;
        } // delta_t is already set
      } else if ( ARG_IS_S_EXACT(0,"compress1") ) {
        if ( ARG_IS_S_EXACT(1,"linear") ) { 
	  compressA_name=strdup(argv[1]);
	}
        else if (ARG_IS_S_EXACT(1,"discard1")) { 
	  compressA_name=strdup(argv[1]);
	}
        else if (ARG_IS_S_EXACT(1,"discard2")) { 
	  compressA_name=strdup(argv[1]);
	}
	else {
	  Tcl_AppendResult(interp, "Compression function ", argv[1], (char *)NULL);
	  Tcl_AppendResult(interp, " is not implemented. ", (char *)NULL);
	  return TCL_ERROR;
	}
        argc -= 2;
        argv += 2; 
      } else if ( ARG_IS_S_EXACT(0,"compress2") ) {
          if ( ARG_IS_S_EXACT(1,"linear") )  { 
	    compressB_name=strdup(argv[1]);
	  }
          else if (ARG_IS_S_EXACT(1,"discard1")) { 
	    compressB_name=strdup(argv[1]);
	  }
          else if (ARG_IS_S_EXACT(1,"discard2")) { 
	    compressB_name=strdup(argv[1]);
	  } else { 
	    Tcl_AppendResult(interp, "Compression function ", argv[1], (char *)NULL); 
	    Tcl_AppendResult(interp, " is not implemented. ", (char *)NULL); 
	    return TCL_ERROR; 
	  }
          argc -= 2;
          argv += 2; 
      } else if ( ARG0_IS_S("update") ) {
        sprintf(buffer,"%d ",no);
        Tcl_AppendResult(interp, "Correlation error: cannot update correlation ", buffer, (char *)NULL);
        Tcl_AppendResult(interp, " It must be set up first", (char *)NULL);
        return TCL_ERROR;
      } else if ( ARG0_IS_S("print") ) {
        sprintf(buffer,"%d ",no);
        Tcl_AppendResult(interp, "Correlation error: cannot print correlation ", buffer, (char *)NULL);
        Tcl_AppendResult(interp, " It must be set up first", (char *)NULL);
        return TCL_ERROR;
      } else {
        Tcl_AppendResult(interp, "unknown argument ", argv[0], (char *)NULL);
        return TCL_ERROR;
      }
    }
  } else {
      sprintf(buffer,"%d \n",no);
      Tcl_AppendResult(interp,"Correlation error: Correlations must be numbered consecutively",buffer,(char *)NULL);
      return TCL_ERROR;
  }


  // Obtain references to observables
  Observable* A=observables.at(obs_A_id);
  Observable* B=observables.at(obs_B_id);
  // Dimensionallity of observables


  // Instantiate and store the correlation object
  correlations.insert(std::make_pair(no,new Correlation(delta_t, tau_lin, tau_max, 1, 
  dim_A, dim_B, dim_corr, *A, *B, 
  corr_operation_name, compressA_name, compressB_name, *args)));


  return TCL_OK;
}


int parse_corr_operation(Tcl_Interp* interp, int argc, char** argv, int* change, char **corr_operation_name, unsigned int* dim_corr, unsigned int dim_A, unsigned int dim_B, void **args) {
  *corr_operation_name = strdup(argv[0]);
  DoubleList *dl;
  int i;
  if (ARG_IS_S_EXACT(0,"componentwise_product")) {
    *change=1;
    return TCL_OK;
  } else if (ARG_IS_S_EXACT(0,"complex_conjugate_product")) {
    *change=1;
    return TCL_OK;
  } else if (ARG_IS_S_EXACT(0,"tensor_product")) {
    *change=1;
    return TCL_OK;
  } else if (ARG_IS_S_EXACT(0,"square_distance_componentwise")) {
    *change=1;
    return TCL_OK;
  } else if (ARG_IS_S_EXACT(0,"fcs_acf")) {
    dl = (DoubleList*)Utils::malloc(sizeof(DoubleList));
    alloc_doublelist(dl, 3);
    dl->n = 3;
    for (i=0; i<3; i++ )  { 
        if ( ! ARG_IS_D(i+1,dl->e[i])  )  { 
            Tcl_AppendResult(interp, "\nError: fcs_acf mut be folowed by 3 arguments of type double, got: ", argv[i+1], "\n\n" , (char *)NULL);
            return TCL_ERROR;
        }
    }
    // convert beam waist radii to their squares
    for (i=0; i<3; i++ ) 
        dl->e[i] *=  dl->e[i];
    *args = (void*) dl;
    *change=4;
    return TCL_OK;
  } else if (ARG0_IS_S("scalar_product")) {
    *change=1;
    return TCL_OK;
  } else {
    Tcl_AppendResult(interp, "Unknown correlation operation: ", argv[0], "\n" , (char *)NULL);
    return TCL_ERROR;
  }
}



int double_correlation_print_correlation( Correlation& corr, Tcl_Interp* interp) {
  double dt=corr.dt;
  char buffer[TCL_DOUBLE_SPACE];
//  char ibuffer[TCL_INTEGER_SPACE+2];

  for (unsigned j=0; j<corr.n_result; j++) {
     Tcl_AppendResult(interp, " { ", (char *)NULL);
     Tcl_PrintDouble(interp, corr.tau[j]*dt, buffer);
     Tcl_AppendResult(interp, buffer, " ",(char *)NULL);
     sprintf(buffer, "%d ", corr.n_sweeps[j]);
     Tcl_AppendResult(interp, buffer, " ",(char *)NULL);
     for (unsigned k=0; k< corr.dim_corr; k++) {
     if (corr.n_sweeps[j] == 0 ) {
       Tcl_PrintDouble(interp, 0., buffer);
       Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
     }
     else {
       Tcl_PrintDouble(interp, corr.result[j][k]/ (double) corr.n_sweeps[j], buffer);
       Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
     }
     }
     Tcl_AppendResult(interp, " } \n", (char *)NULL);
  }
  return 0;
}


int double_correlation_print_spherically_averaged_sf(Correlation& corr, Tcl_Interp* interp) {

  int qi,qj,qk,qn, dim_sf, order2;
  double dt=corr.dt;
  observable_sf_params* params=(observable_sf_params*)corr.A_obs.container;
  char buffer[TCL_DOUBLE_SPACE];
  int *q_vals;
  double *q_density;
  double *av_sf_Re;
  double *av_sf_Im;
  
  q_vals=params->q_vals;
  q_density=params->q_density;
  dim_sf=params->dim_sf;
  order2=params->order*params->order;
  
  av_sf_Re=(double*)Utils::malloc(order2*sizeof(double));
  av_sf_Im=(double*)Utils::malloc(order2*sizeof(double));

  // compute spherically averaged sf
  for (unsigned j=0; j<corr.n_result; j++) {
    // compute the spherically averaged sf for current dt
    for (int k=0;k<order2;k++) av_sf_Re[k]=av_sf_Im[k]=0.0;
    for (int k=0;k<dim_sf;k++) {
      qi=q_vals[3*k  ];
      qj=q_vals[3*k+1];
      qk=q_vals[3*k+2];
      qn= qi*qi + qj*qj + qk*qk;
      av_sf_Re[qn-1]+=corr.result[j][2*k  ];
      av_sf_Im[qn-1]+=corr.result[j][2*k+1];
    }
    for (int k=0;k<order2;k++) { 
      if (q_density[k]>0.0) {
        av_sf_Re[k]/=q_density[k];
        av_sf_Im[k]/=q_density[k];
      }
      // note: if q_density[k]==0, we did not add anything to av_sf_Xx[k], so it is 0.0
    }
    // now print what we obtained
    for (int k=0;k<order2;k++) { 
      Tcl_AppendResult(interp, " { ", (char *)NULL);
      Tcl_PrintDouble(interp, corr.tau[j]*dt, buffer);
      Tcl_AppendResult(interp, buffer, " } { ",(char *)NULL);
      Tcl_PrintDouble(interp, (double) corr.n_sweeps[j], buffer);
      Tcl_AppendResult(interp, buffer, " } { ", (char *)NULL);
      Tcl_PrintDouble(interp, sqrt((double)k+1.0), buffer);
      Tcl_AppendResult(interp, buffer, " } { ", (char *)NULL);
      if (corr.n_sweeps[j] == 0 ) {
        Tcl_PrintDouble(interp, 0., buffer);
        Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
        Tcl_PrintDouble(interp, 0., buffer);
        Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
      }
      else {
        Tcl_PrintDouble(interp, av_sf_Re[k]/(double) corr.n_sweeps[j], buffer);
        Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
        Tcl_PrintDouble(interp, av_sf_Im[k]/(double) corr.n_sweeps[j], buffer);
        Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
      }
      Tcl_AppendResult(interp, " } \n", (char *)NULL);
    }
  }
  return 0;
}

int double_correlation_write_to_file( Correlation& corr, char* filename) {
  FILE* file=0;
  double dt=corr.dt;
  file=fopen(filename, "w");
  if (!file) {
    return 1;
  }
  for (unsigned j=0; j<corr.n_result; j++) {
    fprintf(file, "%.6g %d ", corr.tau[j]*dt, corr.n_sweeps[j]);
    for (unsigned k=0; k< corr.dim_corr; k++) {
      if (corr.n_sweeps[j] == 0 )
        fprintf(file, "%.6g ", 0.);
      else 
        fprintf(file, "%.6g ", corr.result[j][k]/ (double) corr.n_sweeps[j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
  return 0;
}


