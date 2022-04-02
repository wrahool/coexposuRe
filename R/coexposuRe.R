#' Get Simulated Network
#'
#' `simulate_single_network()` runs a single simulation and returns coexposure networks and additional information related to those networks. \cr
#' @references Mukerjee, S. (2021). A systematic comparison of community detection algorithms for measuring selective exposure in co-exposure networks. Scientific Reports, 11(1), 1-11.
#' @param n1 The number of media outlets in the environment
#' @param n2 The number of agents in the environment
#' @param n3 The number of types of media outlets and agents in the environment. Must be less than n1 and n2.
#' @param rho The randomizing parameter. When rho is 0, agents strictly visit outlets of the same type as themselves. When rho is 1, agents visit outlets randomly.
#' @param a the exponent of the power-law distribution from which the reputation of the media outlets is drawn. Default is 2.
#' @param b the skewness parameter that determines how skewed the browsing behavior of the agents is. Can be negative or positive, but not 0. Default is 1.
#' @param show_network if TRUE, will plot the coexposure network. Default is FALSE.
#' @param niter used if `show_network` = TRUE. Passed to `plot.igraph` for visualizing network.
#' @return `simulate_single_network` returns a list with three elements: \cr
#' * `g` the coexposure network of media outlets \cr
#' * `ag` the augmented coexposure network \crggplot2::
#' * `outlet_dat` a tibble containing the details of the simulated media outlets
#' @examples
#' res <- simulate_single_network(n1 = 50, n2 = 30, n3 = 3, rho = 0.1, show_network = TRUE)
#' @export
simulate_single_network <- function(n1, n2, n3, rho, a = 2, b = 1, show_network = FALSE, niter = 500) {

  # check if parameters are valid


  if(n3 > n1) stop("n3 (number of types of outlets/agents) cannot be more than n1 (number of outlets)")
  if(n3 > n2) stop("n3 (number of types of outlets/agents) cannot be more than n2 (number of agents)")
  if(b == 0) stop("b (skewness) cannot be 0")

  if(class(n1) != "numeric" | class(n2) != "numeric" | class(n3) != "numeric" | class(rho) != "numeric") {
    stop("n1, n2, and n3 should be of type numeric.")
  }

  if(class(a) != "numeric" | class(b) != "numeric") {
    stop("a and b should be of type numeric.")
  }

  # stopifnot(!missing(n1), !missing(n2), !missing(n3), !missing(rho))
  # stopifnot(class(n1) == "numeric")
  # stopifnot(class(n2) == "numeric")
  # stopifnot(class(n3) == "numeric")
  # stopifnot(n3 <= n1)
  # stopifnot(n3 <= n2)
  # stopifnot(b != 0)


  message("Initializing universe...")
  outlet_ids <- 1:n1
  p_ids <- 1:n2
  types <- as.character(1:n3)

  if(a > 1) {
    outlet_rep <-  poweRlaw::rpldis(n1, 1, alpha = a) # power law distribution
    outlet_rep_normalized <- outlet_rep / sum(outlet_rep)
  } else {
    outlet_rep <- rep(1, n1)
    outlet_rep_normalized <- outlet_rep / sum(outlet_rep)
  }

  outlets_tbl <- tibble::tibble(
    outlet_id = outlet_ids,
    outlet_name = paste("O", outlet_ids, sep = "_"),
    outlet_type = sample_atleast_once(types, n1), # at least one website of each type
    outlet_repute = outlet_rep_normalized
  )

  all_n4 <- fGarch::rsnorm(n = n2, mean = 0, sd = 1, xi = b)
  all_n4_scaled <- round(((all_n4 - min(all_n4))/(max(all_n4)-min(all_n4))*(n1-1) + 1))

  audience_tbl <- tibble::tibble(
    p_id = p_ids,
    p_name = paste("P", p_ids, sep = ""),
    p_type = sample_atleast_once(types, n2),       # at least one audience member of each type
    p_n4 = all_n4_scaled
  )

  audience_el <- NULL

  # loop over each person
  message("Simulating audience behavior...")
  for(p in 1:n2) {

    n4 <- audience_tbl$p_n4[p]

    # when rho is 0 all of their choices are selective
    # when rho is 1 all of their choices are random
    random_choices_allowed <- round(rho * n4)
    selective_choices_allowed <- n4 - random_choices_allowed

    selective_outlets_pool <- outlets_tbl %>%
      dplyr::filter(outlet_type == audience_tbl$p_type[p]) %>%
      dplyr::select(outlet_id, outlet_repute)

    if(nrow(selective_outlets_pool) == 1) {                         # this is to prevent the bug where 1 type gets 1 outlet
      selective_outlets_pool <- selective_outlets_pool %>%
        rbind(selective_outlets_pool)
    }

    selective_chosen_outlets <- selective_outlets_pool %>%          # from
      dplyr::pull(outlet_id) %>%                                           # all outlets in the selective outlets pool
      sample(selective_choices_allowed, replace = TRUE,             # randomly sample with prob = outlet repute (R auto-normalizes the probabilities of the subset)
             prob = selective_outlets_pool$outlet_repute)

    random_chosen_outlets <- outlets_tbl %>%                        # from
      dplyr::pull(outlet_id) %>%                                           # all outlets in the universe
      sample(random_choices_allowed, replace = TRUE,                # randomly sample with prob = outlet_repute
             prob = outlets_tbl$outlet_repute)

    all_chosen_outlets <- c(selective_chosen_outlets,
                            random_chosen_outlets)

    # build the edge-list for the audience network
    audience_el <- audience_el %>%
      rbind(
        tibble::tibble(
          p_name = paste("P", rep(p, n4), sep = "_"),               # one column is the p_id
          outlet_name = paste("O", all_chosen_outlets, sep = "_")   # second column is the outlet_id
        )
      )
  }

  outlet_reach <- audience_el %>%
    dplyr::pull(outlet_name) %>%
    table() %>%
    tibble::as_tibble() %>%
    dplyr::rename(uv = n) %>%
    dplyr::select(outlet_name = 1, everything())

  message("Constructing network...")
  audience_g <- igraph::graph_from_data_frame(audience_el, directed = F)
  igraph::V(audience_g)$type <- substr(igraph::V(audience_g)$name, 1, 1) == "O"

  # graph projections
  projection_graphs <- igraph::bipartite_projection(audience_g, multiplicity = TRUE)
  outlet_projection <- projection_graphs$proj2

  igraph::V(outlet_projection)$type <- igraph::V(outlet_projection)$name %>%
    lapply(FUN = function(x) {
      outlets_tbl %>%
        dplyr::filter(outlet_name == x) %>%
        dplyr::pull(outlet_type)
    }
    ) %>%
    unlist()

  # # colrs <- c("gray50", "tomato", "gold", "purple", "cyan")
  # V(outlet_projection)$color <- ifelse(V(outlet_projection)$type == "A", "gray50",
  #                                      ifelse(V(outlet_projection)$type == "B", "tomato",
  #                                             ifelse(V(outlet_projection)$type == "C", "gold",
  #                                                    ifelse(V(outlet_projection)$type == "D", "olivedrab4",
  #                                                           "cyan"))))

  outlet_projection_sl <- outlet_projection
  outlet_projection_sl[from=igraph::V(outlet_projection_sl), to=igraph::V(outlet_projection_sl)] = 1
  for(v in igraph::V(outlet_projection_sl)$name) {
    igraph::E(outlet_projection_sl)[v %--% v]$weight <- outlet_reach %>%
      dplyr::filter(outlet_name == v) %>%
      dplyr::pull(uv)
  }

  if(show_network) {

    l <- igraph::layout_with_fr(outlet_projection, niter=niter)
    if(n3 < 10) {
      pal <- RColorBrewer::brewer.pal(n3, "Set1")
    } else {
      pal <- randomcoloR::randomColor(n3)
    }

    igraph::plot.igraph(outlet_projection,
                        edge.color="gray90",
                        vertex.color = pal[as.numeric(as.factor(igraph::vertex_attr(outlet_projection, "type")))],
                        vertex.size = 5,
                        layout = l,
                        vertex.label = NA)

  }

  return(list("g" = outlet_projection, "ag" = outlet_projection_sl, "outlet_data" = outlets_tbl))
}

#' Generate and analyze simulated networks
#'
#' Calls `simulate_single_network` with different values of rho and returns the results of the analysis of those selected networks
#' Two analyses are currently supported with support for more planned for the future. See the `analyze` parameter for details.
#' @references Mukerjee, S. (2021). A systematic comparison of community detection algorithms for measuring selective exposure in co-exposure networks. Scientific Reports, 11(1), 1-11.
#' @param n1 the number of media outlets in the environment
#' @param n2 the number of agents in the environment
#' @param n3 the number of types of media outlets and agents in the environment. Must be less than n1 and n2.
#' @param rho_min the minimum value of the randomizing parameter to be used in the simulation. When rho is 0, agents strictly visit outlets of the same type as themselves. When rho is 1, agents visit outlets randomly.
#' @param rho_max the maximum value of the randomizing parameter to be used in the simulation. When rho is 0, agents strictly visit outlets of the same type as themselves. When rho is 1, agents visit outlets randomly.
#' @param rho_inc the value by which rho changes for each step of the simulation
#' @param a the exponent of the power-law distribution from which the reputation of the media outlets is drawn. Default is 2.
#' @param b the skewness parameter that determines how skewed the browsing behavior of the agents is. Can be negative or positive, but not 0. Default is 1.
#' @param N number of networks to build for each value of rho. Defaults to 100.
#' @param analyze what analysis to perform? Currently supports two: "cd" for community detection. "c" for centralization.
#' @param metric only applicable if analyze == "cd". Which variant of NMI to use for calculation. Default is "max". Can be "max", "min", "sqrt", "sum", "joint"
#' @param correct only applicable if analyze == "cd". If TRUE, NMI scores are scaled using the scaling factor. This corrects for overfitting by community detection algorithms. Default is TRUE.
#' @param plot_results if TRUE, a graph of the results is shown.
#' @return a tibble with the results for each iteration within the simulation
#' @examples
#' res <- analyze_simulated_networks(n1 = 50, n2 = 30, n3 = 3, rho_min = 0, rho_max = 0.4, N = 5, analyze = "c" )
#' res <- analyze_simulated_networks(n1 = 50, n2 = 30, n3 = 3, rho_min = 0, rho_max = 0.4, N = 5, analyze = "cd" , correct = F)
#' @export
analyze_simulated_networks <- function(n1, n2, n3, rho_min = 0, rho_max = 1, rho_inc = 0.1, a = 2, b = 1, N = 100, analyze, metric = "max", correct = T, plot_results = T) {

  if (analyze == "cd") {
    res_tbl <- NULL

    for(rho in seq(from = rho_min, to = rho_max, by = rho_inc)) {
      i <- 1
      while(i <= N) {
        message(paste0("rho : ", rho, " iteration : ", i))

        sim_res <- simulate_single_network(n1, n2, n3, rho = rho, a, b)

        g <- sim_res$g
        g_sl <- sim_res$ag
        o_tbl <- sim_res$outlet_data

        if(length(igraph::V(g)) <= 1) {
          # message("Rerun...")
          next
        }

        c_wt <- tryCatch(
          # cluster_walktrap uses E(g)$weight by default
          igraph::cluster_walktrap(g),
          error = function(e) {
            return(NA)
          })

        c_wt2 <- tryCatch(
          igraph::cluster_walktrap(g_sl),
          error = function(e) {
            return(NA)
          })

        c_l <- tryCatch(
          # cluster_louvain uses E(g)$weight by default
          igraph::cluster_louvain(g),
          error = function(e) {
            return(NA)
          })

        c_l2 <- tryCatch(
          # cluster_louvain uses E(g)$weight by default
          igraph::cluster_louvain(g_sl),
          error = function(e) {
            return(NA)
          })

        c_fg <- tryCatch(
          # cluster_fast_greedy uses E(g)$weight by default
          igraph::cluster_fast_greedy(g),
          error = function(e) {
            return(NA)
          })

        c_fg2 <- tryCatch(
          # cluster_fast_greedy uses E(g)$weight by default
          igraph::cluster_fast_greedy(g_sl),
          error = function(e) {
            return(NA)
          })

        c_eb <- tryCatch(
          # cluster edge_betweenness uses E(g)$weight by default
          igraph::cluster_edge_betweenness(g),
          error = function(e) {
            return(NA)
          })

        c_eb2 <- tryCatch(
          # cluster edge_betweenness uses E(g)$weight by default
          igraph::cluster_edge_betweenness(g_sl),
          error = function(e) {
            return(NA)
          })

        c_im <- tryCatch(
          # cluster_infomap needs an argument called e.weights, but uses E(g)$weight by default
          igraph::cluster_infomap(g),
          error = function(e) {
            return(NA)
          })

        c_im2 <- tryCatch(
          # cluster_infomap needs an argument called e.weights, but uses E(g)$weight by default
          igraph::cluster_infomap(g_sl),
          error = function(e) {
            return(NA)
          })


        c_lp <- tryCatch(
          # label propagation uses weight by default
          igraph::cluster_label_prop(g),
          error = function(e) {
            return(NA)
          })

        c_lp2 <- tryCatch(
          # label propagation uses weight by default
          igraph::cluster_label_prop(g_sl),
          error = function(e) {
            return(NA)
          })

        c_le <- tryCatch(
          # leading eigenvector uses weight by default
          igraph::cluster_leading_eigen(g, options = list(maxiter=1000000)),
          error = function(e) {
            return(NA)
          })

        c_le2 <- tryCatch(
          # leading eigenvector uses E(g)$weight by default
          igraph::cluster_leading_eigen(g_sl, options = list(maxiter=1000000)),
          error = function(e) {
            return(NA)
          })

        c_sl <- tryCatch(
          # for spinglass, by default weights = NULL and that uses the E(g)$weight attribute
          igraph::cluster_spinglass(g),
          error = function(e) {
            return(NA)
          })

        c_sl2 <- tryCatch(
          # for spinglass, by default weights = NULL and that uses the E(g)$weight attribute
          igraph::cluster_spinglass(g_sl),
          error = function(e) {
            return(NA)
          })


        all_cs <- list(c_wt, c_wt2,
                       c_l, c_l2,
                       c_fg, c_fg2,
                       c_eb, c_eb2,
                       c_im, c_im2,
                       c_lp, c_lp2,
                       c_le, c_le2,
                       c_sl, c_sl2
        )

        cd_used <- c(
          "wt", "wt2",
          "l", "l2",
          "fg", "fg2",
          "eb", "eb2",
          "im", "im2",
          "lp", "lp2",
          "le", "le2",
          "sl", "sl2"
        )

        NMI_scores <- sapply(all_cs, FUN = function(x) {
          get_NMI(x, o_tbl, metric)
        })

        if(correct) {
          scaling_factors <- sapply(all_cs, FUN = function(x) {
            get_scalingfactor(x, o_tbl)
          })

          res_tbl <- tibble::tibble(
            run = i,
            rho = rho,
            method = cd_used,
            NMI_scores = NMI_scores,
            scaling_factors = scaling_factors
          ) %>%
            rbind(res_tbl)

        } else {

          res_tbl <- tibble::tibble(
            run = i,
            rho = rho,
            method = cd_used,
            NMI_scores = NMI_scores
          ) %>%
            rbind(res_tbl)
        }

        i <- i+1
      }
    }

    res_tbl <- res_tbl %>%
      dplyr::mutate(network_type = ifelse(grepl("2", method), "augmented", "baseline")) %>%
      dplyr::mutate(method = gsub("2", "", method)) %>%
      dplyr::mutate(method = ifelse(method == "wt", "Walktrap",
                             ifelse(method == "l", "Multilevel",
                                    ifelse(method == "fg", "Fast Greedy",
                                           ifelse(method == "eb", "Edge Betweenness",
                                                  ifelse(method == "im", "Infomap",
                                                         ifelse(method == "lp", "Label Propagation",
                                                                ifelse(method == "le", "Leading Eigenvector",
                                                                       ifelse(method == "sl", "Spin Glass", "Unknown")))))))))


    if(correct) {

      message("Scaling NMI values...")
      res_tbl <- res_tbl %>%
        dplyr::mutate(SNMI_scores = NMI_scores * scaling_factors)

      suppressMessages(plot_tbl <- res_tbl %>%
        dplyr::group_by(method, network_type, rho) %>%
        dplyr::summarise(meanSNMI = mean(SNMI_scores),
                  sdSNMI = sd(SNMI_scores)))

      if(plot_results) {

        message("Plotting results...")
        print(ggplot2::ggplot(data = plot_tbl, ggplot2::aes(x = rho, fill = network_type, color = network_type)) +
                ggplot2::geom_line(ggplot2::aes(y = meanSNMI)) +
                ggplot2::labs(y = "scaled NMI") +
                ggplot2::geom_ribbon(ggplot2::aes(ymin = meanSNMI - sdSNMI, ymax = meanSNMI + sdSNMI), alpha = 0.3) +
                ggplot2::facet_wrap(~method, nrow = 2))
      }

    } else {

      message("NMI values not scaled.")

      suppressMessages(plot_tbl <- res_tbl %>%
        dplyr::group_by(method, network_type, rho) %>%
        dplyr::summarise(meanNMI = mean(NMI_scores),
                  sdNMI = sd(NMI_scores)))

      if(plot_results) {

        message("Plotting results...")
        print(ggplot2::ggplot(data = plot_tbl, ggplot2::aes(x = rho, fill = network_type, color = network_type)) +
                ggplot2::geom_line(ggplot2::aes(y = meanNMI)) +
                ggplot2::geom_ribbon(ggplot2::aes(ymin = meanNMI - sdNMI, ymax = meanNMI + sdNMI), alpha = 0.3) +
                ggplot2::facet_wrap(~method, nrow = 2))
      }

    }

    return(res_tbl)

  }

  if(analyze == "c") {

    res_tbl <- NULL

    for(rho in seq(from = rho_min, to = rho_max, by = rho_inc)) {

      i <- 1

      while(i <= N) {
        message(paste0("rho : ", rho, " iteration : ", i))

        sim_res <- simulate_single_network(n1, n2, n3, rho = rho, a, b)

        g <- sim_res$g
        g_sl <- sim_res$ag
        o_tbl <- sim_res$outlet_data

        if(length(igraph::V(g)) <= 1) {
          # message("Rerun...")
          next
        }

        c_d <- tryCatch(
          igraph::centr_degree(g, normalized = F)$centralization,
          error = function(e) {
            NA
          })

        c_d_n <- tryCatch(
          igraph::centr_degree(g, normalized = T)$centralization,
          error = function(e) {
            NA
          })

        c_b <- tryCatch(
          igraph::centr_betw(g, normalized = F)$centralization,
          error = function(e) {
            NA
          })


        c_b_n <- tryCatch(
          igraph::centr_betw(g, normalized = T)$centralization,
          error = function(e) {
            NA
          })

        c_c <- tryCatch(
          igraph::centr_clo(g, normalized = F)$centralization,
          error = function(e) {
            NA
          })

        c_c_n <- tryCatch(
          igraph::centr_clo(g, normalized = T)$centralization,
          error = function(e) {
            NA
          })

        c_e <- tryCatch(
          igraph::centr_eigen(g, normalized = F)$centralization,
          error = function(e) {
            NA
          })

        c_e_n <- tryCatch(
          igraph::centr_eigen(g, normalized = T)$centralization,
          error = function(e) {
            NA
          })

        res_tbl <- tibble::tibble(
          run = i,
          rho = rho,
          c_deg = c_d,
          c_deg_norm = c_d_n,
          c_bet = c_b,
          c_bet_norm = c_b_n,
          c_clos = c_c,
          c_clos_norm = c_c_n,
          c_eigen = c_e,
          c_eigen_norm = c_e_n
        ) %>%
          rbind(res_tbl)

        i <- i+1
      }
    }

    res_tbl <- res_tbl %>%
      tidyr::pivot_longer(cols = starts_with("c_"), names_to = "type", values_to = "centralization")

    if(plot_results) {
      plot_tbl <- res_tbl %>%
        dplyr::group_by(rho, type) %>%
        dplyr::summarize(mean_c = mean(centralization, na.rm = T),
                  sd_c = sd(centralization, na.rm = T)) %>%
        dplyr::mutate(type = ifelse(type == "c_deg", "Regular degree",
                             ifelse(type == "c_deg_norm", "Normalized degree",
                                    ifelse(type == "c_bet", "Regular betweenness",
                                           ifelse(type == "c_bet_norm", "Normalized betweenness",
                                                  ifelse(type == "c_clos", "Regular closeness",
                                                         ifelse(type == "c_clos_norm", "Normalized closeness",
                                                                ifelse(type == "c_eigen", "Regular eigenvector", "Normalized eigenvector")))))))) %>%
        tidyr::separate(type, into = c("type1", "type2"), sep = " ")

      print(ggplot2::ggplot(data = plot_tbl, ggplot2::aes(x = rho)) +
              ggplot2::geom_line(ggplot2::aes(y = mean_c)) +
              ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_c - sd_c, ymax = mean_c + sd_c), alpha = 0.3) +
              ggplot2::labs(y = "centralization") +
              ggplot2::facet_wrap(type1 ~ type2, scales = "free", nrow = 2))
    }

    return(res_tbl)
  }
}

