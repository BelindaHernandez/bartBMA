# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

find_internal_nodes_pred <- function(treetable) {
    .Call(`_bartBMA_find_internal_nodes_pred`, treetable)
}

find_term_nodes_pred <- function(tree_table) {
    .Call(`_bartBMA_find_term_nodes_pred`, tree_table)
}

get_original_pred <- function(low, high, sp_low, sp_high, sum_preds) {
    .Call(`_bartBMA_get_original_pred`, low, high, sp_low, sp_high, sum_preds)
}

bartBMA_get_testdata_term_obs_pred <- function(test_data, tree_data, term_node_means) {
    .Call(`_bartBMA_bartBMA_get_testdata_term_obs_pred`, test_data, tree_data, term_node_means)
}

get_BART_BMA_test_predictions <- function(test_data, BIC, sum_trees, y_minmax) {
    .Call(`_bartBMA_get_BART_BMA_test_predictions`, test_data, BIC, sum_trees, y_minmax)
}

get_imp_vars <- function(split_vars, num_col, current_vars) {
    .Call(`_bartBMA_get_imp_vars`, split_vars, num_col, current_vars)
}

get_weighted_var_imp <- function(num_vars, BIC, sum_trees) {
    .Call(`_bartBMA_get_weighted_var_imp`, num_vars, BIC, sum_trees)
}

csample_num <- function(x, size, replace, prob = as.numeric( c())) {
    .Call(`_bartBMA_csample_num`, x, size, replace, prob)
}

add_rows <- function(prior_tree_table_temp, grow_node) {
    .Call(`_bartBMA_add_rows`, prior_tree_table_temp, grow_node)
}

addcol <- function(prior_tree_matrix_temp, grow_node, ld_obs, rd_obs) {
    .Call(`_bartBMA_addcol`, prior_tree_matrix_temp, grow_node, ld_obs, rd_obs)
}

set_daughter_to_end_tree <- function(grow_node, prior_tree_table_temp, left_daughter) {
    .Call(`_bartBMA_set_daughter_to_end_tree`, grow_node, prior_tree_table_temp, left_daughter)
}

set_daughter_to_end_mat <- function(d, prior_tree_matrix_temp, left_daughter, ld_obs, rd_obs) {
    .Call(`_bartBMA_set_daughter_to_end_mat`, d, prior_tree_matrix_temp, left_daughter, ld_obs, rd_obs)
}

remove_zero <- function(nodes_at_depth) {
    .Call(`_bartBMA_remove_zero`, nodes_at_depth)
}

order_intvec_ <- function(x) {
    .Call(`_bartBMA_order_intvec_`, x)
}

get_gnp <- function(nodes_at_depth, grow_node) {
    .Call(`_bartBMA_get_gnp`, nodes_at_depth, grow_node)
}

find_term_nodes <- function(tree_table) {
    .Call(`_bartBMA_find_term_nodes`, tree_table)
}

find_term_obs <- function(tree_matrix_temp, terminal_node) {
    .Call(`_bartBMA_find_term_obs`, tree_matrix_temp, terminal_node)
}

likelihood_function <- function(y_temp, treetable_temp, obs_to_nodes_temp, a, mu, nu, lambda) {
    .Call(`_bartBMA_likelihood_function`, y_temp, treetable_temp, obs_to_nodes_temp, a, mu, nu, lambda)
}

find_internal_nodes <- function(treetable) {
    .Call(`_bartBMA_find_internal_nodes`, treetable)
}

find_prev_nonterm <- function(find_nonterm, prev) {
    .Call(`_bartBMA_find_prev_nonterm`, find_nonterm, prev)
}

find_nodes_to_update <- function(all_ld, left_daughter) {
    .Call(`_bartBMA_find_nodes_to_update`, all_ld, left_daughter)
}

set_tree_to_middle <- function(node_to_update, prior_tree_table_temp, grow_node, left_daughter) {
    .Call(`_bartBMA_set_tree_to_middle`, node_to_update, prior_tree_table_temp, grow_node, left_daughter)
}

update_grow_obs <- function(prior_tree_matrix_temp, grow_node, left_daughter, d, ld_obs, rd_obs) {
    .Call(`_bartBMA_update_grow_obs`, prior_tree_matrix_temp, grow_node, left_daughter, d, ld_obs, rd_obs)
}

find_obs_to_update_grow <- function(prior_tree_matrix_temp, left_daughter, d, ld_obs, rd_obs) {
    .Call(`_bartBMA_find_obs_to_update_grow`, prior_tree_matrix_temp, left_daughter, d, ld_obs, rd_obs)
}

get_subset <- function(xmat, grow_obs) {
    .Call(`_bartBMA_get_subset`, xmat, grow_obs)
}

get_daughter_obs <- function(xmat, obs_to_update, split_var, split_point) {
    .Call(`_bartBMA_get_daughter_obs`, xmat, obs_to_update, split_var, split_point)
}

find_term_cols <- function(tree_matrix_temp, terminal_node) {
    .Call(`_bartBMA_find_term_cols`, tree_matrix_temp, terminal_node)
}

get_grow_obs <- function(xmat, grow_obs, split_var) {
    .Call(`_bartBMA_get_grow_obs`, xmat, grow_obs, split_var)
}

grow_tree <- function(xmat, y, prior_tree_matrix, grow_node, prior_tree_table, splitvar, splitpoint, terminal_nodes, grow_obs, d, get_min, data_curr_node) {
    .Call(`_bartBMA_grow_tree`, xmat, y, prior_tree_matrix, grow_node, prior_tree_table, splitvar, splitpoint, terminal_nodes, grow_obs, d, get_min, data_curr_node)
}

set_daughter <- function(left_daughter, right_daughter, ld_obs, rd_obs, tree_matrix_temp, term_cols) {
    .Call(`_bartBMA_set_daughter`, left_daughter, right_daughter, ld_obs, rd_obs, tree_matrix_temp, term_cols)
}

order_ <- function(x) {
    .Call(`_bartBMA_order_`, x)
}

get_tree_prior <- function(tree_table, tree_matrix, alpha, beta) {
    .Call(`_bartBMA_get_tree_prior`, tree_table, tree_matrix, alpha, beta)
}

start_tree <- function(start_mean, start_sd) {
    .Call(`_bartBMA_start_tree`, start_mean, start_sd)
}

start_matrix <- function(n) {
    .Call(`_bartBMA_start_matrix`, n)
}

evaluate_model_occams_window <- function(tree_lik, lowest_BIC, c, tree_list, tree_mat_list, tree_parent) {
    .Call(`_bartBMA_evaluate_model_occams_window`, tree_lik, lowest_BIC, c, tree_list, tree_mat_list, tree_parent)
}

get_testdata_term_obs <- function(test_data, tree_data, term_node_means) {
    .Call(`_bartBMA_get_testdata_term_obs`, test_data, tree_data, term_node_means)
}

resize <- function(x, n) {
    .Call(`_bartBMA_resize`, x, n)
}

resize_bigger <- function(x, n) {
    .Call(`_bartBMA_resize_bigger`, x, n)
}

J <- function(treetable_temp, obs_to_nodes_temp, tree_term_nodes) {
    .Call(`_bartBMA_J`, treetable_temp, obs_to_nodes_temp, tree_term_nodes)
}

mu_vector <- function(sum_treetable, n) {
    .Call(`_bartBMA_mu_vector`, sum_treetable, n)
}

W <- function(sum_treetable, sum_obs_to_nodes, n) {
    .Call(`_bartBMA_W`, sum_treetable, sum_obs_to_nodes, n)
}

sumtree_likelihood_function <- function(y_temp, sum_treetable, sum_obs_to_nodes, n, a, nu, lambda) {
    .Call(`_bartBMA_sumtree_likelihood_function`, y_temp, sum_treetable, sum_obs_to_nodes, n, a, nu, lambda)
}

get_best_split <- function(resids, data, treetable, tree_mat, a, mu, nu, lambda, c, lowest_BIC, parent, cp_mat, alpha, beta, maxOWsize, first_round) {
    .Call(`_bartBMA_get_best_split`, resids, data, treetable, tree_mat, a, mu, nu, lambda, c, lowest_BIC, parent, cp_mat, alpha, beta, maxOWsize, first_round)
}

get_best_split_sum <- function(resids, data, treetable, tree_mat, a, mu, nu, lambda, c, lowest_BIC, parent, cp_mat, alpha, beta, maxOWsize, first_round, sum_trees, sum_trees_mat, y_scaled, parent2, i) {
    .Call(`_bartBMA_get_best_split_sum`, resids, data, treetable, tree_mat, a, mu, nu, lambda, c, lowest_BIC, parent, cp_mat, alpha, beta, maxOWsize, first_round, sum_trees, sum_trees_mat, y_scaled, parent2, i)
}

update_mean_var <- function(tree_table, tree_matrix, resids, a) {
    .Call(`_bartBMA_update_mean_var`, tree_table, tree_matrix, resids, a)
}

update_predictions <- function(tree_table, tree_matrix, new_mean, n) {
    .Call(`_bartBMA_update_predictions`, tree_table, tree_matrix, new_mean, n)
}

subsetter <- function(a, b) {
    .Call(`_bartBMA_subsetter`, a, b)
}

order_inc_ <- function(x) {
    .Call(`_bartBMA_order_inc_`, x)
}

min_which2 <- function(array, n, minout, whichout) {
    .Call(`_bartBMA_min_which2`, array, n, minout, whichout)
}

mll_meanvar2 <- function(x, x2, n) {
    .Call(`_bartBMA_mll_meanvar2`, x, x2, n)
}

PELT_meanvar_norm2 <- function(resp, pen) {
    .Call(`_bartBMA_PELT_meanvar_norm2`, resp, pen)
}

SS <- function(x, y, split) {
    .Call(`_bartBMA_SS`, x, y, split)
}

gridCP <- function(x, y, gridSize = 10L) {
    .Call(`_bartBMA_gridCP`, x, y, gridSize)
}

make_gridpoint_cpmat <- function(data, resp, gridsize, num_cp) {
    .Call(`_bartBMA_make_gridpoint_cpmat`, data, resp, gridsize, num_cp)
}

make_pelt_cpmat <- function(data, resp, pen, num_cp) {
    .Call(`_bartBMA_make_pelt_cpmat`, data, resp, pen, num_cp)
}

get_best_trees <- function(D1, resids, a, mu, nu, lambda, c, sigma_mu, tree_table, tree_mat, lowest_BIC, first_round, parent, cp_mat_list, err_list, test_data, alpha, beta, is_test_data, pen, num_cp, split_rule_node, gridpoint, maxOWsize) {
    .Call(`_bartBMA_get_best_trees`, D1, resids, a, mu, nu, lambda, c, sigma_mu, tree_table, tree_mat, lowest_BIC, first_round, parent, cp_mat_list, err_list, test_data, alpha, beta, is_test_data, pen, num_cp, split_rule_node, gridpoint, maxOWsize)
}

get_best_trees_sum <- function(D1, resids, a, mu, nu, lambda, c, sigma_mu, tree_table, tree_mat, lowest_BIC, first_round, parent, cp_mat_list, err_list, test_data, alpha, beta, is_test_data, pen, num_cp, split_rule_node, gridpoint, maxOWsize, prev_sum_trees, prev_sum_trees_mat, y_scaled) {
    .Call(`_bartBMA_get_best_trees_sum`, D1, resids, a, mu, nu, lambda, c, sigma_mu, tree_table, tree_mat, lowest_BIC, first_round, parent, cp_mat_list, err_list, test_data, alpha, beta, is_test_data, pen, num_cp, split_rule_node, gridpoint, maxOWsize, prev_sum_trees, prev_sum_trees_mat, y_scaled)
}

scale_response <- function(a, b, c, d, y) {
    .Call(`_bartBMA_scale_response`, a, b, c, d, y)
}

get_original <- function(low, high, sp_low, sp_high, sum_preds) {
    .Call(`_bartBMA_get_original`, low, high, sp_low, sp_high, sum_preds)
}

BART_BMA_sumLikelihood <- function(data, y, start_mean, start_sd, a, mu, nu, lambda, c, sigma_mu, pen, num_cp, test_data, num_rounds, alpha, beta, split_rule_node, gridpoint, maxOWsize) {
    .Call(`_bartBMA_BART_BMA_sumLikelihood`, data, y, start_mean, start_sd, a, mu, nu, lambda, c, sigma_mu, pen, num_cp, test_data, num_rounds, alpha, beta, split_rule_node, gridpoint, maxOWsize)
}

find_term_nodes_gs <- function(tree_table) {
    .Call(`_bartBMA_find_term_nodes_gs`, tree_table)
}

find_term_obs_gs <- function(tree_matrix_temp, terminal_node) {
    .Call(`_bartBMA_find_term_obs_gs`, tree_matrix_temp, terminal_node)
}

calc_rowsums <- function(predictions) {
    .Call(`_bartBMA_calc_rowsums`, predictions)
}

calculate_resids <- function(predictions, response) {
    .Call(`_bartBMA_calculate_resids`, predictions, response)
}

update_Gibbs_mean_var <- function(tree_table, resids, a, sigma, mu_mu, terminal_nodes, term_obs_tree) {
    .Call(`_bartBMA_update_Gibbs_mean_var`, tree_table, resids, a, sigma, mu_mu, terminal_nodes, term_obs_tree)
}

update_sigma <- function(a1, b, resids, n) {
    .Call(`_bartBMA_update_sigma`, a1, b, resids, n)
}

find_node_means <- function(sum_tree, term_nodes) {
    .Call(`_bartBMA_find_node_means`, sum_tree, term_nodes)
}

get_tree_info <- function(overall_sum_trees, overall_sum_mat, num_obs) {
    .Call(`_bartBMA_get_tree_info`, overall_sum_trees, overall_sum_mat, num_obs)
}

remove_curr_col <- function(predy, i) {
    .Call(`_bartBMA_remove_curr_col`, predy, i)
}

get_new_mean <- function(terminal_nodes, new_mean_var) {
    .Call(`_bartBMA_get_new_mean`, terminal_nodes, new_mean_var)
}

update_predictions_gs <- function(tree_table, new_mean, new_var, n, terminal_nodes, term_obs_tree) {
    .Call(`_bartBMA_update_predictions_gs`, tree_table, new_mean, new_var, n, terminal_nodes, term_obs_tree)
}

scale_response_gs <- function(a, b, c, d, y) {
    .Call(`_bartBMA_scale_response_gs`, a, b, c, d, y)
}

get_original_gs <- function(low, high, sp_low, sp_high, sum_preds) {
    .Call(`_bartBMA_get_original_gs`, low, high, sp_low, sp_high, sum_preds)
}

find_internal_nodes_gs <- function(treetable) {
    .Call(`_bartBMA_find_internal_nodes_gs`, treetable)
}

get_tree_info_test_data <- function(test_data, tree_data) {
    .Call(`_bartBMA_get_tree_info_test_data`, test_data, tree_data)
}

get_tree_info_testdata_overall <- function(overall_sum_trees, num_obs, test_data) {
    .Call(`_bartBMA_get_tree_info_testdata_overall`, overall_sum_trees, num_obs, test_data)
}

gibbs_sampler <- function(overall_sum_trees, overall_sum_mat, y, BIC_weights, num_iter, burnin, num_obs, num_test_obs, a, sigma, mu_mu, nu, lambda, resids, test_data) {
    .Call(`_bartBMA_gibbs_sampler`, overall_sum_trees, overall_sum_mat, y, BIC_weights, num_iter, burnin, num_obs, num_test_obs, a, sigma, mu_mu, nu, lambda, resids, test_data)
}

gibbs_sampler2 <- function(overall_sum_trees, overall_sum_mat, y, BIC_weights, num_iter, burnin, num_obs, a, sigma, mu_mu, nu, lambda, resids) {
    .Call(`_bartBMA_gibbs_sampler2`, overall_sum_trees, overall_sum_mat, y, BIC_weights, num_iter, burnin, num_obs, a, sigma, mu_mu, nu, lambda, resids)
}

