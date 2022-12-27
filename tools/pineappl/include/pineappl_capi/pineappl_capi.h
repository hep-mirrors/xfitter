/*
 * PineAPPL - PDF-independent binning of phase space weights
 * Copyright (C) 2020-2021  Christopher Schwan

 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PINEAPPL_H
#define PINEAPPL_H

/* Generated with cbindgen:0.24.3 */


#define PINEAPPL_CAPI_MAJOR 0
#define PINEAPPL_CAPI_MINOR 5
#define PINEAPPL_CAPI_PATCH 8


#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/**
 * Main data structure of `PineAPPL`. This structure contains a `Subgrid` for each `LumiEntry`,
 * bin, and coupling order it was created with.
 */
typedef struct pineappl_grid pineappl_grid;

/**
 * Key-value storage for passing optional information during grid creation with
 * `pineappl_grid_new`.
 */
typedef struct pineappl_keyval pineappl_keyval;

/**
 * Type for defining a luminosity function.
 */
typedef struct pineappl_lumi pineappl_lumi;

/**
 * Type for reading and accessing subgrids.
 */
typedef struct pineappl_subgrid pineappl_subgrid;

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/**
 * Returns the number of bins in `grid`.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call.
 */
size_t pineappl_grid_bin_count(const struct pineappl_grid *grid);

/**
 * Returns the number of dimensions of the bins this grid has.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call.
 */
size_t pineappl_grid_bin_dimensions(const struct pineappl_grid *grid);

/**
 * Stores the bin sizes of `grid` in `bin_sizes`.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call. The parameter `bin_sizes` must point to an array that is as
 * long as `grid` has bins.
 */
void pineappl_grid_bin_normalizations(const struct pineappl_grid *grid, double *bin_sizes);

/**
 * Write the left limits for the specified dimension into `left`.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call. The parameter `left` must point to an array that is as large
 * as `grid` has bins. If `dimension` is larger or equal the number of dimensions for this grid,
 * nothing is written into `left`, the result is undefined.
 */
void pineappl_grid_bin_limits_left(const struct pineappl_grid *grid,
                                   size_t dimension,
                                   double *left);

/**
 * Write the right limits for the specified dimension into `right`.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call. The parameter `right` must point to an array that is as
 * large as `grid` has bins. If `dimension` is larger or equal the number of dimensions for this
 * grid, nothing is written into `right`, the result is undefined.
 */
void pineappl_grid_bin_limits_right(const struct pineappl_grid *grid,
                                    size_t dimension,
                                    double *right);

/**
 * Returns a cloned object of `grid`.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call.
 */
struct pineappl_grid *pineappl_grid_clone(const struct pineappl_grid *grid);

/**
 * Convolutes the specified grid with the PDF `xfx`, which is the PDF for a hadron with the PDG id
 * `pdg_id`, and strong coupling `alphas`. These functions must evaluate the PDFs for the given
 * `x` and `q2` for the parton with the given PDG id, `pdg_id`, and return the result. Note that
 * the value must be the PDF multiplied with its argument `x`. The value of the pointer `state`
 * provided to these functions is the same one given to this function. The parameter `order_mask`
 * must be as long as there are perturbative orders contained in `grid` and is used to selectively
 * disable (`false`) or enable (`true`) individual orders. If `order_mask` is set to `NULL`, all
 * orders are active. The parameter `lumi_mask` can be used similarly, but must be as long as the
 * luminosity function `grid` was created with has entries, or `NULL` to enable all luminosities.
 * The values `xi_ren` and `xi_fac` can be used to vary the renormalization and factorization from
 * its central value, which corresponds to `1.0`. After convolution of the grid with the PDFs the
 * differential cross section for each bin is written into `results`.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call. The function pointers `xfx1`, `xfx2`, and `alphas` must not
 * be null pointers and point to valid functions. The parameters `order_mask` and `lumi_mask` must
 * either be null pointers or point to arrays that are as long as `grid` has orders and lumi
 * entries, respectively. Finally, `results` must be as long as `grid` has bins.
 */
void pineappl_grid_convolute_with_one(const struct pineappl_grid *grid,
                                      int32_t pdg_id,
                                      double (*xfx)(int32_t pdg_id, double x, double q2, void *state),
                                      double (*alphas)(double q2, void *state),
                                      void *state,
                                      const bool *order_mask,
                                      const bool *lumi_mask,
                                      double xi_ren,
                                      double xi_fac,
                                      double *results);

/**
 * Convolutes the specified grid with the PDFs `xfx1` and `xfx2`, which are the PDFs of hadrons
 * with PDG ids `pdg_id1` and `pdg_id2`, respectively, and strong coupling `alphas`. These
 * functions must evaluate the PDFs for the given `x` and `q2` for the parton with the given PDG
 * id, `pdg_id`, and return the result. Note that the value must be the PDF multiplied with its
 * argument `x`. The value of the pointer `state` provided to these functions is the same one
 * given to this function. The parameter `order_mask` must be as long as there are perturbative
 * orders contained in `grid` and is used to selectively disable (`false`) or enable (`true`)
 * individual orders. If `order_mask` is set to `NULL`, all orders are active. The parameter
 * `lumi_mask` can be used similarly, but must be as long as the luminosity function `grid` was
 * created with has entries, or `NULL` to enable all luminosities. The values `xi_ren` and
 * `xi_fac` can be used to vary the renormalization and factorization from its central value,
 * which corresponds to `1.0`. After convolution of the grid with the PDFs the differential cross
 * section for each bin is written into `results`.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call. The function pointers `xfx1`, `xfx2`, and `alphas` must not
 * be null pointers and point to valid functions. The parameters `order_mask` and `lumi_mask` must
 * either be null pointers or point to arrays that are as long as `grid` has orders and lumi
 * entries, respectively. Finally, `results` must be as long as `grid` has bins.
 */
void pineappl_grid_convolute_with_two(const struct pineappl_grid *grid,
                                      int32_t pdg_id1,
                                      double (*xfx1)(int32_t pdg_id, double x, double q2, void *state),
                                      int32_t pdg_id2,
                                      double (*xfx2)(int32_t pdg_id, double x, double q2, void *state),
                                      double (*alphas)(double q2, void *state),
                                      void *state,
                                      const bool *order_mask,
                                      const bool *lumi_mask,
                                      double xi_ren,
                                      double xi_fac,
                                      double *results);

/**
 * Delete a grid previously created with `pineappl_grid_new`.
 */
void pineappl_grid_delete(struct pineappl_grid *grid);

/**
 * Fill `grid` for the given momentum fractions `x1` and `x2`, at the scale `q2` for the given
 * value of the `order`, `observable`, and `lumi` with `weight`.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call.
 */
void pineappl_grid_fill(struct pineappl_grid *grid,
                        double x1,
                        double x2,
                        double q2,
                        size_t order,
                        double observable,
                        size_t lumi,
                        double weight);

/**
 * Fill `grid` for the given momentum fractions `x1` and `x2`, at the scale `q2` for the given
 * value of the `order` and `observable` with `weights`. The parameter of weight must contain a
 * result for entry of the luminosity function the grid was created with.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call.
 */
void pineappl_grid_fill_all(struct pineappl_grid *grid,
                            double x1,
                            double x2,
                            double q2,
                            size_t order,
                            double observable,
                            const double *weights);

/**
 * Fill `grid` with as many points as indicated by `size`.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call. Additionally, all remaining pointer parameters must be
 * arrays as long as specified by `size`.
 */
void pineappl_grid_fill_array(struct pineappl_grid *grid,
                              const double *x1,
                              const double *x2,
                              const double *q2,
                              const size_t *orders,
                              const double *observables,
                              const size_t *lumis,
                              const double *weights,
                              size_t size);

/**
 * Return the luminosity function of `grid`.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call.
 */
struct pineappl_lumi *pineappl_grid_lumi(const struct pineappl_grid *grid);

/**
 * Write the order parameters of `grid` into `order_params`.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call. The pointer `order_params` must point to an array as large
 * as four times the number of orders in `grid`.
 */
void pineappl_grid_order_params(const struct pineappl_grid *grid, uint32_t *order_params);

/**
 * Return the number of orders in `grid`.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call.
 */
size_t pineappl_grid_order_count(const struct pineappl_grid *grid);

/**
 * Creates a new and empty grid. The creation requires four different sets of parameters:
 * - The luminosity function `lumi`: A pointer to the luminosity function that specifies how the
 * cross section should be reconstructed.
 * - Order specification `orders` and `order_params`. Each `PineAPPL` grid contains a number of
 * different perturbative orders, specified by `orders`. The array `order_params` stores the
 * exponent of each perturbative order and must contain 4 integers denoting the exponent of the
 * string coupling, of the electromagnetic coupling, of the logarithm of the renormalization
 * scale, and finally of the logarithm of the factorization scale.
 * - The observable definition `bins` and `bin_limits`. Each `PineAPPL` grid can store observables
 * from a one-dimensional distribution. To this end `bins` specifies how many observables are
 * stored and `bin_limits` must contain `bins + 1` entries denoting the left and right limit for
 * each bin.
 * - More (optional) information can be given in a key-value storage `key_vals`, which might be
 * a null pointer, to signal there are no further parameters that need to be set.
 *
 * # Safety
 *
 * The parameter `lumi` must point a valid luminosity function created by `pineappl_lumi_new`.
 * `order_params` must be an array with a length of `4 * orders`, and `bin_limits` an array with
 * length `bins + 1`. `key_vals` must be a valid `KeyVal` object created by `pineappl_keyval_new`.
 *
 * # Panics
 *
 * TODO
 */
struct pineappl_grid *pineappl_grid_new(const struct pineappl_lumi *lumi,
                                        size_t orders,
                                        const uint32_t *order_params,
                                        size_t bins,
                                        const double *bin_limits,
                                        const struct pineappl_keyval *key_vals);

/**
 * Read a `PineAPPL` grid from a file with name `filename`.
 *
 * # Safety
 *
 * The parameter `filename` must be a C string pointing to an existing grid file.
 *
 * # Panics
 *
 * TODO
 */
struct pineappl_grid *pineappl_grid_read(const char *filename);

/**
 * Merges `other` into `grid` and subsequently deletes `other`.
 *
 * # Safety
 *
 * Both `grid` and `other` must be valid `Grid` objects created by either `pineappl_grid_new` or
 * `pineappl_grid_read`. If `other` is a `NULL` pointer, this function does not do anything.
 *
 * # Panics
 *
 * TODO
 */
void pineappl_grid_merge_and_delete(struct pineappl_grid *grid, struct pineappl_grid *other);

/**
 * Scale all grids in `grid` by `factor`.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call.
 */
void pineappl_grid_scale(struct pineappl_grid *grid, double factor);

/**
 * Optimizes the grid representation for space efficiency.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call.
 */
void pineappl_grid_optimize(struct pineappl_grid *grid);

/**
 * Scales each subgrid by a factor which is the product of the given values `alphas`, `alpha`,
 * `logxir`, and `logxif`, each raised to the corresponding powers for each subgrid. In addition,
 * every subgrid is scaled by a factor `global` independently of its order.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call.
 */
void pineappl_grid_scale_by_order(struct pineappl_grid *grid,
                                  double alphas,
                                  double alpha,
                                  double logxir,
                                  double logxif,
                                  double global);

/**
 * Return the value for `key` stored in `grid`. If `key` isn't found, `NULL` will be returned.
 * After usage the string must be deallocated using [`pineappl_string_delete`].
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the `NULL`
 * pointer, this function is not safe to call. The parameter `key` must be non-`NULL` and a valid
 * C string.
 *
 * # Panics
 *
 * TODO
 */
char *pineappl_grid_key_value(const struct pineappl_grid *grid, const char *key);

/**
 * Sets an internal key-value pair for the grid.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call. The parameters `key` and `value` must be non-`NULL` and
 * valid C strings.
 *
 * # Panics
 *
 * TODO
 */
void pineappl_grid_set_key_value(struct pineappl_grid *grid, const char *key, const char *value);

/**
 * Sets a remapper for the grid. This can be used to 'upgrade' one-dimensional bin limits to
 * N-dimensional ones. The new bin limits must be given in the form of tuples giving the left and
 * right limits, and a tuple for each dimension.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call. The arrays `normalizations` and `limits` must be at least as
 * long as the number of bins of the grid and `2 * dimensions * bins`, respectively.
 *
 * # Panics
 *
 * TODO
 */
void pineappl_grid_set_remapper(struct pineappl_grid *grid,
                                size_t dimensions,
                                const double *normalizations,
                                const double *limits);

/**
 * Write `grid` to a file with name `filename`. If `filename` ends in `.lz4` the grid is
 * automatically LZ4 compressed.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call. The parameter `filename` must be a non-`NULL`, non-empty,
 * and valid C string pointing to a non-existing, but writable file.
 *
 * # Panics
 *
 * TODO
 */
void pineappl_grid_write(const struct pineappl_grid *grid, const char *filename);

/**
 * Adds a linear combination of initial states to the luminosity function `lumi`.
 *
 * # Safety
 *
 * The parameter `lumi` must point to a valid `Lumi` object created by `pineappl_lumi_new`.
 * `pdg_id_pairs` must be an array with length `2 * combinations`, and `factors` with length of
 * `combinations`.
 */
void pineappl_lumi_add(struct pineappl_lumi *lumi,
                       size_t combinations,
                       const int32_t *pdg_id_pairs,
                       const double *factors);

/**
 * Returns the number of combinations of the luminosity function `lumi` for the specified entry.
 *
 * # Safety
 *
 * The parameter `lumi` must point to a valid `Lumi` object created by `pineappl_lumi_new` or
 * `pineappl_grid_lumi`.
 */
size_t pineappl_lumi_combinations(const struct pineappl_lumi *lumi, size_t entry);

/**
 * Returns the number of channel for the luminosity function `lumi`.
 *
 * # Safety
 *
 * The parameter `lumi` must point to a valid `Lumi` object created by `pineappl_lumi_new` or
 * `pineappl_grid_lumi`.
 */
size_t pineappl_lumi_count(const struct pineappl_lumi *lumi);

/**
 * Delete luminosity function previously created with `pineappl_lumi_new`.
 */
void pineappl_lumi_delete(struct pineappl_lumi *lumi);

/**
 * Read out the channel with index `entry` of the luminosity function `lumi`. The PDG ids and
 * factors will be copied into `pdg_ids` and `factors`.
 *
 * # Safety
 *
 * The parameter `lumi` must point to a valid `Lumi` object created by `pineappl_lumi_new` or
 * `pineappl_grid_lumi`. The parameter `factors` must point to an array as long as the size
 * returned by `pineappl_lumi_combinations` and `pdg_ids` must point to an array that is twice as
 * long.
 */
void pineappl_lumi_entry(const struct pineappl_lumi *lumi,
                         size_t entry,
                         int32_t *pdg_ids,
                         double *factors);

/**
 * Creates a new luminosity function and returns a pointer to it. If no longer needed, the object
 * should be deleted using `pineappl_lumi_delete`.
 */
struct pineappl_lumi *pineappl_lumi_new(void);

/**
 * Fills `buffer` with the q2-slice for the index `q2_slice` of the grid for the specified
 * `order`, `bin`, and `lumi`.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call. `buffer` must be as large as the square of the return value
 * of `pineappl_subgrid_x_grid_count`.
 */
void pineappl_grid_export_mu2_slice(const struct pineappl_grid *grid,
                                    size_t order,
                                    size_t bin,
                                    size_t lumi,
                                    size_t q2_slice,
                                    double *buffer);

/**
 * Write into `tuple` the lower and upper limit of filled q2 slices for the grid with the
 * specified indices. The last slice that is filled is one minus the upper limit.
 *
 * # Safety
 *
 * If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
 * this function is not safe to call. `tuple` must point to an array with two elements.
 */
void pineappl_grid_nonzero_mu2_slices(const struct pineappl_grid *grid,
                                      size_t order,
                                      size_t bin,
                                      size_t lumi,
                                      size_t *tuple);

/**
 * Deletes a subgrid created with [`pineappl_subgrid_new2`]. If `subgrid` is the null pointer,
 * nothing is done.
 */
void pineappl_subgrid_delete(struct pineappl_subgrid *subgrid);

/**
 * This function takes replaces the subgrid in `grid` with the corresponding indices `order`,
 * `bin` and `lumi` with the one given in `subgrid`. If `subgrid` is the null pointer, the specied
 * subgrid is replaced with an empty one.
 *
 * # Safety
 *
 * Both `grid` and `subgrid` must point to valid objects. The parameter `subgrid` can be the null
 * pointer.
 */
void pineappl_grid_replace_and_delete(struct pineappl_grid *grid,
                                      struct pineappl_subgrid *subgrid,
                                      size_t order,
                                      size_t bin,
                                      size_t lumi);

/**
 * Creates a new subgrid, which can be filled with [`pineappl_subgrid_import_mu2_slice`]. The
 * array `mu2_grid` must contain the (squared) values of the renormalization and then the
 * factorization scale, such that twice the value of `mu2_grid_len` gives the length of the array.
 *
 * # Safety
 *
 * The pointers `mu2_grid`, `x1_grid` and `x2_grid` must be non-`NULL` and array. Furthermore, the
 * array `mu2_grid` must be *twice* as long as given in `mu2_grid_len`, and the arrays `x1_grid`
 * and `x2_grid` as long as specified by `x1_grid_len` and `x2_grid_len`, respectively.
 */
struct pineappl_subgrid *pineappl_subgrid_new2(size_t mu2_grid_len,
                                               const double *mu2_grid,
                                               size_t x1_grid_len,
                                               const double *x1_grid,
                                               size_t x2_grid_len,
                                               const double *x2_grid);

/**
 * Imports `slice` for the given index into `subgrid`.
 *
 * # Safety
 *
 * The parameter `subgrid` and the array `slice` must be non-`NULL` and `slice` must be at least
 * as long as the product `x1_grid_len * x2_grid_len` that were used to create the subgrid with.
 * The index `mu2_slice` must be smaller than `mu2_grid_len` that was used in
 * [`pineappl_subgrid_new2`] to create `subgrid`.
 */
void pineappl_subgrid_import_mu2_slice(struct pineappl_subgrid *subgrid,
                                       size_t mu2_slice,
                                       const double *slice);

/**
 * Delete the previously created object pointed to by `key_vals`.
 */
void pineappl_keyval_delete(struct pineappl_keyval *key_vals);

/**
 * Get the boolean-valued parameter with name `key` stored in `key_vals`.
 *
 * # Safety
 *
 * The parameter `key_vals` must point to a valid `KeyVal` object created by
 * `pineappl_keyval_new`. `key` must be a valid C string.
 */
bool pineappl_keyval_bool(const struct pineappl_keyval *key_vals, const char *key);

/**
 * Get the double-valued parameter with name `key` stored in `key_vals`.
 *
 * # Safety
 *
 * The parameter `key_vals` must point to a valid `KeyVal` object created by
 * `pineappl_keyval_new`. `key` must be a valid C string.
 */
double pineappl_keyval_double(const struct pineappl_keyval *key_vals, const char *key);

/**
 * Get the string-valued parameter with name `key` stored in `key_vals`.
 *
 * # Safety
 *
 * The parameter `key_vals` must point to a valid `KeyVal` object created by
 * `pineappl_keyval_new`. `key` must be a valid C string.
 */
int32_t pineappl_keyval_int(const struct pineappl_keyval *key_vals, const char *key);

/**
 * Get the int-valued parameter with name `key` stored in `key_vals`.
 *
 * # Safety
 *
 * The parameter `key_vals` must point to a valid `KeyVal` object created by
 * `pineappl_keyval_new`. `key` must be a valid C string.
 */
const char *pineappl_keyval_string(const struct pineappl_keyval *key_vals, const char *key);

/**
 * Return a pointer to newly-created `pineappl_keyval` object.
 */
struct pineappl_keyval *pineappl_keyval_new(void);

/**
 * Set the double-valued parameter with name `key` to `value` in `key_vals`.
 *
 * # Safety
 *
 * The parameter `key_vals` must point to a valid `KeyVal` object created by
 * `pineappl_keyval_new`. `key` must be a valid C string.
 */
void pineappl_keyval_set_bool(struct pineappl_keyval *key_vals, const char *key, bool value);

/**
 * Set the double-valued parameter with name `key` to `value` in `key_vals`.
 *
 * # Safety
 *
 * The parameter `key_vals` must point to a valid `KeyVal` object created by
 * `pineappl_keyval_new`. `key` must be a valid C string.
 */
void pineappl_keyval_set_double(struct pineappl_keyval *key_vals, const char *key, double value);

/**
 * Set the int-valued parameter with name `key` to `value` in `key_vals`.
 *
 * # Safety
 *
 * The parameter `key_vals` must point to a valid `KeyVal` object created by
 * `pineappl_keyval_new`. `key` must be a valid C string.
 */
void pineappl_keyval_set_int(struct pineappl_keyval *key_vals, const char *key, int32_t value);

/**
 * Set the string-valued parameter with name `key` to `value` in `key_vals`.
 *
 * # Safety
 *
 * The parameter `key_vals` must point to a valid `KeyVal` object created by
 * `pineappl_keyval_new`. `key` must be a valid C string.
 */
void pineappl_keyval_set_string(struct pineappl_keyval *key_vals,
                                const char *key,
                                const char *value);

/**
 * Deletes a string previously allocated by [`pineappl_grid_key_value`]. If `string` is a `NULL`
 * pointer this function does nothing.
 *
 * # Safety
 *
 * The parameter `string` must be a pointer to string created by [`pineappl_grid_key_value`] or
 * `NULL`, otherwise this function is not safe to call.
 */
void pineappl_string_delete(char *string);

#ifdef __cplusplus
} // extern "C"
#endif // __cplusplus

#endif /* PINEAPPL_H */
