#include <Rcpp.h>
/*
#include "minc2.h"
#include "minc_cpp.h"*/

#include "minc2-simple.h"

#include <sstream>
#include <stdlib.h>
using namespace Rcpp;
using namespace std;


static vector<minc2_file_handle> _open_minc2_volumes(CharacterVector filenames){
  
  vector<minc2_file_handle> volumes;
  CharacterVector::iterator file_iterator;
  
  for(file_iterator = filenames.begin();
      file_iterator != filenames.end();
      ++file_iterator)
  {

      minc2_file_handle h=minc2_allocate0();
      if(minc2_open(h,*file_iterator)!=MINC2_SUCCESS)
          throw "Can't open file";
      volumes.push_back(h);
  }
  
  return(volumes);
}

static void _close_minc2_volumes(vector<minc2_file_handle>& volumes)
{
  vector<minc2_file_handle>::iterator it;
  for(it = volumes.begin();
      it != volumes.end();
      ++it)
      {
            minc2_close(*it);
            minc2_free(*it);
      }
}

// [[Rcpp::export]]
List rcpp_minc_apply2(CharacterVector filenames,
                     bool use_mask,
                     CharacterVector mask,
                     double mask_lower_val,
                     double mask_upper_val,
                     RObject value_for_mask,
                     bool filter_masked,
                     NumericVector slab_sizes,
                     Function fun, 
                     List args) {

  /*using idea from https://github.com/vfonov/minc2-simple/blob/master/example/RInside/example_minc_rinside.cpp*/
  std::vector<minc2_file_handle> volumes;
  minc2_file_handle mask_handle=NULL;
  minc2_file_iterator_handle input_minc_it=NULL;
  minc2_file_iterator_handle mask_minc_it=NULL;
  size_t added_results = 0;
  Rcpp::List results;
  NumericVector result_inds;
  
  try {
    volumes =_open_minc2_volumes(filenames);
    std::vector<double> voxels(volumes.size());  
    
    if(use_mask)
      mask_handle = _open_minc2_volumes(mask)[0]; /*TODO: check volume dimensions*/
    
  
    int nvols = volumes.size();
    int nvoxels = 0;
    minc2_nelement(volumes[0],&nvoxels);
    
    results=Rcpp::List(nvoxels);
    result_inds=Rcpp::NumericVector(nvoxels);
    // Setup looping constructs
    size_t voxel_pos = 0;
    
  
    Rprintf("Number of Volumes: %d Number of Voxels: %d \n", nvols,nvoxels);
    
    input_minc_it=minc2_iterator_allocate0();
    
    
    if(minc2_multi_iterator_input_start(input_minc_it,&volumes[0],MINC2_DOUBLE,volumes.size())!=MINC2_SUCCESS)
      throw "Error in iterators"; /*this will also check that all dimension are the same*/
    
    if(use_mask) {
      mask_minc_it=minc2_iterator_allocate0();
      if(minc2_iterator_input_start(mask_minc_it,mask_handle,MINC2_DOUBLE)!=MINC2_SUCCESS)
        throw "Error in mask iterators"; /*this will also check that all dimension are the same*/
    }
    bool advance=false;
    do {
       bool process_voxel=1;
       //Check for user breaking the loop
       checkUserInterrupt();
       
       //Check if a mask was supplied, then check if the current voxel is masked
       if(use_mask){
         double mask_val;
         minc2_iterator_get_values(mask_minc_it,&mask_val); /*TODO: check return code*/
          
         process_voxel=( (mask_val > (mask_lower_val - .5) && mask_val < (mask_upper_val + .5))); /*why .5 ?*/
       }
       
       if(process_voxel) {
         minc2_iterator_get_values(input_minc_it,&voxels[0]); /*TODO: check return code*/
        results[added_results] = fun(voxels, args);
        result_inds[added_results] = voxel_pos;
        ++added_results;
       }
       
       voxel_pos++; /*TODO: change to the index pos?*/
       if(use_mask)
         minc2_iterator_next(mask_minc_it);
    } while(minc2_iterator_next(input_minc_it)==MINC2_SUCCESS);
  
  
    Rprintf("\n");
  } catch(const char *err) {
    if(use_mask && mask_minc_it!=NULL)
      minc2_iterator_free(mask_minc_it);
    if(input_minc_it!=NULL)
      minc2_iterator_free(input_minc_it);
    _close_minc2_volumes(volumes);
    stop(err);
  }

  if(use_mask && mask_minc_it!=NULL)
    minc2_iterator_free(mask_minc_it);
  if(input_minc_it!=NULL)
    minc2_iterator_free(input_minc_it);
  _close_minc2_volumes(volumes);
  
  results.erase(results.begin() + added_results, results.end());
  result_inds.erase(result_inds.begin() + added_results, result_inds.end());
  
  return List::create(_["vals"] = results,
                      _["inds"] = result_inds);
}

/* kate: indent-mode cstyle; indent-width 2; replace-tabs on; */
