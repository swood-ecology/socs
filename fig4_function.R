#### Effect plot ##############################################################
#### For generating Figure 4
#### Modified from jtools package

#' @param robust Should robust standard errors be used to find confidence
#'   intervals for supported models? Default is FALSE, but you should specify
#'   the type of sandwich standard errors if you'd like to use them (i.e.,
#'   `"HC0"`, `"HC1"`, and so on). If `TRUE`, defaults to `"HC3"` standard
#'   errors.
#'
#' @param cluster For clustered standard errors, provide the column name of
#'   the cluster variable in the input data frame (as a string). Alternately,
#'   provide a vector of clusters.
#'
#' @param pred.values Values of `pred` to use instead of the equi-spaced
#'   series by default (for continuous variables) or all unique values (for
#'   non-continuous variables).
#' @param force.cat Force the continuous `pred` to be treated as categorical?
#'   default is FALSE, but this can be useful for things like dummy 0/1
#'   variables.

effect_plot <- function(model, pred, pred.values = NULL, data = NULL,
                        robust = FALSE, cluster = NULL, x.label = NULL,
                        y.label = NULL, colors = NULL, shape = NULL,
                        force.cat = FALSE, ...) {

  # Evaluate the pred arg
  pred <- quo_name(enexpr(pred))

  # Have a sensible interval default for categorical predictors
  if ("interval" %nin% names(match.call())[-1] &
    !(is.numeric(get_data(model, warn = FALSE)[[pred]]) &
      force.cat == FALSE)) {
    interval <- TRUE
  }

  if (force.cat == TRUE & is.null(pred.values)) {
    if (is.null(data)) {
      data <- get_data(model)
    }
    pred.values <- sort(unique(suppressMessages(data[[pred]])))
  }

  # Deal with legacy color argument
  pred_out <- jtools::make_predictions(model,
    pred = pred, pred.values = pred.values,
    at = NULL, center = "all",
    interval = TRUE, int.type = c("confidence", "prediction"),
    outcome.scale = "response", robust = robust,
    cluster = cluster, vcov = NULL,
    set.offset = 1, return.orig.data = TRUE,
    data = data, ...
  )

  # Putting these outputs into separate objects
  pm <- pred_out[[1]]
  d <- pred_out[[2]]

  if (is.numeric(d[[pred]]) & force.cat == FALSE) {
    plot_effect_continuous(
      predictions = pm, pred = pred,
      data = d, x.label = x.label, y.label = y.label,
      resp = get_response_name(model),
      colors = colors, shape = shape,
      weights = get_weights(model, d)$weights_name
    )
  }
}

plot_effect_continuous <-
  function(predictions, pred,
             data = NULL, x.label = NULL, y.label = NULL,
             colors = NULL, shape = NULL, resp = NULL, weights = NULL) {
    pm <- predictions
    d <- data

    if (is.null(x.label)) {
      x.label <- pred
    }

    if (is.null(y.label)) {
      y.label <- resp
    }

    pred <- sym(pred)
    resp <- sym(resp)
    if (!is.null(weights)) {
      weights <- sym(weights)
    }

    # Starting plot object
    p <- ggplot(pm, aes(x = !!pred, y = !!resp))

    # Define line thickness
    p <- p + geom_path(size = 1, alpha = 0.5)

    # Plot ribbon around line
    p <- p + geom_ribbon(
      data = pm,
      aes(ymin = !!sym("ymin"), ymax = !!sym("ymax")),
      alpha = 1 / 10, show.legend = TRUE
    )

    # Plot points
    p <- p + geom_point(
      data = d,
      aes(
        x = !!pred, y = !!resp, fill = !!colors,
        shape = !!shape
      ),
      size = 2
    )

    # Change color and shape aesthetics
    p <- p + scale_shape_manual(
      name = "Replicate",
      values = c(21, 22, 23, 24)
    ) +
      scale_fill_manual(
        name = "System",
        values = c(
          "#d73027", "#ffff33", "#a6d96a",
          "#006837", "#2166ac"
        )
      ) +
        guides(fill = guide_legend(override.aes = list(shape = 21)))

    # Plot rug
    p <- p + geom_rug(
      data = d,
      mapping = aes(x = !!pred, y = !!resp), alpha = 0.6,
      sides = "lb", inherit.aes = TRUE
    )

    # Using theme_apa for theming...but using legend title and side positioning
    p <- p + theme_minimal()

    p <- p + labs(x = x.label, y = y.label) # better labels for axes

    # Return the plot
    return(p)
  }
