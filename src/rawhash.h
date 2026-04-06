#ifndef RAWHASH_H
#define RAWHASH_H

#include "rutils.h"

/*************
 * options   *
 *************/

#include "roptions.h"

/**************
 * Signal	  *
 * Processing *
 **************/

#include "rsig.h"

/*************
 * index     *
 *************/

#include "rindex.h"

/*************
 * sketch     *
 *************/

#include "rsketch.h"

/*************
 * mapping   *
 *************/

#include "ketopt.h"
#include "rmap.h"

struct config_t {
 ketopt_t o = KETOPT_INIT;
 ri_mapopt_t opt = {0};
 ri_idxopt_t ipt = {0};
 char *fnw = nullptr, *fpore = nullptr;
 int n_threads = 3, io_n_threads = 1;
 bool valid = false;
};

#endif // RAWHASH_H
