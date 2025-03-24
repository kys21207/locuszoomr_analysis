#' Locus scatter plot with credible set highlighting
#' 
#' @author kys21207
#' @date 2025-03-22 17:02:23 UTC
#'
#' Produces a base graphics scatter plot from a 'locus' class object. This
#' function is called by [locus_plot()] to generate the scatter plot portion.
#' Can be used manually with [set_layers()].
#'
#' @param loc Object of class 'locus' to use for plot. See [locus].
#' @param index_snp Specifies index SNP or vector of SNPs to show in different colour/symbol.
#' @param pcutoff Cut-off for p value significance (default: 5e-08).
#' @param scheme Vector of 3 colours: normal points, significant points, index SNP(s).
#' @param cs_color Color for credible set variants (default: 'orange').
#' @param cex Base point size.
#' @param cs_cex Size multiplier for credible set variants (default: 1.5).
#' @param cex.axis Font size for axis numbering.
#' @param cex.lab Font size for axis titles.
#' @param xlab x axis title.
#' @param ylab y axis title.
#' @param ylim y axis limits (y1, y2).
#' @param ylim2 Secondary y axis limits for recombination line.
#' @param yzero Force y axis limit to include y=0.
#' @param xticks Show x axis numbers and title.
#' @param border Show bounding box.
#' @param showLD Show LD with colours.
#' @param LD_scheme Vector of colours for LD plotting.
#' @param recomb_col Colour for recombination rate line.
#' @param recomb_offset Offset for recombination line.
#' @param legend_pos Legend position.
#' @param labels SNPs to label.
#' @param label_x,label_y Label positions.
#' @param eqtl_gene Column for eQTL genes.
#' @param beta Column for beta coefficient direction.
#' @param add Add to existing plot.
#' @param align Set par() alignment.
#' @param ... Other plot() arguments.
#' @return No return value. Produces a scatter plot.
#' @importFrom graphics par legend mtext axis points text lines
#' @export

scatter_plot_cs <- function(loc,
                        index_snp = loc$index_snp,
                        pcutoff = 5e-08,
                        scheme = c('grey', 'dodgerblue', 'red'),
                        cs_color = 'orange',
                        cex = 1,
                        cs_cex = 1.5,
                        cex.axis = 0.9,
                        cex.lab = 1,
                        xlab = NULL,
                        ylab = NULL,
                        ylim = NULL,
                        ylim2 = c(0, 100),
                        yzero = (loc$yvar == "logP"),
                        xticks = TRUE,
                        border = FALSE,
                        showLD = TRUE,
                        LD_scheme = c('grey', 'royalblue', 'cyan2', 'green3', 
                                    'orange', 'red', 'purple'),
                        recomb_col = "blue",
                        recomb_offset = 0,
                        legend_pos = 'topleft',
                        labels = NULL,
                        label_x = 4,
                        label_y = 4,
                        eqtl_gene = NULL,
                        beta = NULL,
                        add = FALSE,
                        align = TRUE,
                        ...) {
    
    # Input validation
    if (!inherits(loc, "locus")) stop("Object of class 'locus' required")
    if (is.null(loc$data)) stop("No data points, only gene tracks")
    
    # Store call and initialize data
    .call <- match.call()
    data <- loc$data
    
    # Set axis labels
    if (is.null(xlab)) xlab <- paste("Chromosome", loc$seqname, "(Mb)")
    if (is.null(ylab)) {
        ylab <- if (loc$yvar == "logP") expression(-log[10] ~ "P") else loc$yvar
    }
    
    # Check for LD data
    hasLD <- "ld" %in% colnames(data)
    
    # Set background colors
    if (!("bg" %in% colnames(data))) {
        if (showLD & hasLD) {
            data$bg <- cut(data$ld, -1:6/5, labels = FALSE)
            data$bg[is.na(data$bg)] <- 1L
            data$bg[data[, loc$labs] %in% index_snp] <- 7L
            data <- data[order(data$bg), ]
            LD_scheme <- rep_len(LD_scheme, 7)
            data$bg <- LD_scheme[data$bg]
        } else if (!is.null(eqtl_gene)) {
            bg <- data[, eqtl_gene]
            bg[data[, loc$p] > pcutoff] <- "ns"
            bg <- relevel(factor(bg, levels = unique(bg)), "ns")
            if (is.null(.call$scheme)) scheme <- eqtl_scheme(nlevels(bg))
            data$bg <- scheme[bg]
        } else {
            # Default coloring including CS variants
            data$bg <- scheme[1]  # Default color
            if (loc$yvar == "logP") {
                data$bg[data[, loc$p] < pcutoff] <- scheme[2]
            }
            data$bg[data[, loc$labs] %in% index_snp] <- scheme[3]
            
            # Color CS variants in orange
            if ("cs_num" %in% colnames(data)) {
                data$bg[!is.na(data$cs_num)] <- cs_color
            }
        }
    }
    
    # Set point shapes and sizes
    pch <- rep(21L, nrow(data))
    point_cex <- rep(cex, nrow(data))
    
    # Handle credible sets
    has_cs <- "cs_num" %in% colnames(data)
    if (has_cs) {
        cs_numbers <- unique(data$cs_num[!is.na(data$cs_num)])
        cs_shapes <- c(22L, 24L, 25L, 23L) # square, up tri, down tri, diamond

    # Initialize cs_genename vector
    cs_genename <- character(length(cs_numbers))
      
        for (i in seq_along(cs_numbers)) {
            cs_idx <- which(data$cs_num == cs_numbers[i])
            pch[cs_idx] <- cs_shapes[(i-1) %% length(cs_shapes) + 1]
            point_cex[cs_idx] <- cex * cs_cex
            cs_genename[i] <- unique(data$flames[cs_idx])
        }
    }
   
    # Handle index SNPs and beta coefficients
    pch[data[, loc$labs] %in% index_snp] <- 23L
    if (!is.null(beta)) {
        sig <- data[, loc$p] < pcutoff
        pch[sig] <- 24 + (1 - sign(data[sig, beta])) / 2
    }
    
    # Override with custom attributes if present
    if ("pch" %in% colnames(data)) pch <- data$pch
    col <- if ("col" %in% colnames(data)) data$col else "black"
    if ("cex" %in% colnames(data)) point_cex <- data$cex
    
    # Set plot parameters
    recomb <- !is.null(loc$recomb) & !is.na(recomb_col)
    if (align) {
        op <- par(mar = c(ifelse(xticks, 3, 0.1), 3.5, 2,
                         ifelse(recomb, 3.5, 1.5)))
        on.exit(par(op))
    }
    
    # Set y-axis limits
    if (is.null(ylim)) {
        ylim <- range(data[, loc$yvar], na.rm = TRUE)
        if (yzero & is.null(ylim)) ylim[1] <- min(c(0, ylim[1]))
        if (!is.null(labels) & (border | recomb)) {
            ylim[2] <- ylim[2] + diff(ylim) * 0.08
        }
    }
    
    # Handle recombination offset
    yd <- diff(ylim)
    if (recomb && recomb_offset != 0) ylim[1] <- ylim[1] - yd * recomb_offset
    
    # Define panel.first expression
    panel.first <- quote({
        if (loc$yvar == "logP" & !is.null(pcutoff)) {
            abline(h = -log10(pcutoff), col = 'darkgrey', lty = 2)
        }
        if (recomb) {
            yd2 <- diff(ylim2)
            fy2 <- function(yy) (yy - ylim2[1]) / yd2 * yd + ylim[1]
            ry <- fy2(loc$recomb$value)
            lines(loc$recomb$start, ry, col = recomb_col)
            labs2 <- pretty(ylim2)
            at <- fy2(labs2)
            axis(4, at = at, labels = labs2,
                 las = 1, tcl = -0.3, mgp = c(1.7, 0.5, 0),
                 cex.axis = cex.axis)
            mtext("Recombination rate (%)", 4, cex = cex.lab * par("cex"),
                  line = 1.7,
                  adj = max(c(0.5 - recomb_offset / 2, 0)))
        }
    })
    
    # Add points to existing plot if requested
    if (add) {
        plot.args <- list(x = data[, loc$pos], y = data[, loc$yvar],
                         pch = pch, bg = data$bg, cex = point_cex)
        if (length(new.args <- list(...))) {
            plot.args[names(new.args)] <- new.args
        }
        return(do.call("points", plot.args))
    }
    
    # Create new plot
    bty <- if (border | recomb) 'o' else 'l'
    plot.args <- list(
        x = data[, loc$pos], y = data[, loc$yvar],
        pch = pch, bg = data$bg, col = col,
        las = 1, font.main = 1,
        cex = point_cex, cex.axis = cex.axis, cex.lab = cex.lab,
        xlim = loc$xrange,
        ylim = ylim,
        xlab = if (xticks) xlab else "",
        ylab = ylab,
        bty = bty,
        xaxt = 'n',
        tcl = -0.3,
        mgp = c(1.7, 0.5, 0),
        panel.first = panel.first
    )
    
    if (recomb && recomb_offset != 0) {
        plot.args$ylab <- ""
        plot.args$yaxt <- 'n'
    }
    
    if (length(new.args <- list(...))) {
        plot.args[names(new.args)] <- new.args
    }
    
    do.call("plot", plot.args)
    
    # Add offset y1 axis ticks if needed
    if (recomb && recomb_offset != 0) {
        ypretty <- pretty(c(min(data[, loc$yvar], na.rm = TRUE), ylim[2]))
        axis(2, at = ypretty, las = 1, 
             mgp = c(1.7, 0.5, 0), cex.axis = cex.axis,
             tcl = -0.3)
        mtext(ylab, 2, cex = cex.lab * par("cex"), line = 1.7,
              adj = min(c(0.5 + recomb_offset / 2.7, 1)))
    }
    
    # Add labels if requested
    if (!is.null(labels)) {
        i <- grep("index", labels, ignore.case = TRUE)
        if (length(i) > 0) {
            if (length(index_snp) == 1) {
                labels[i] <- index_snp
            } else {
                labels <- labels[-i]
                labels <- c(index_snp, labels)
            }
        }
        ind <- match(labels, data[, loc$labs])
        if (any(is.na(ind))) {
            message("label ", paste(labels[is.na(ind)], collapse = ", "),
                    " not found")
        }
        lx <- data[ind, loc$pos]
        ly <- data[ind, loc$yvar]
        labs <- data[ind, loc$labs]
        add_labels(lx, ly, labs, label_x, label_y, cex = cex.axis * 0.95)
    }
    
    # Add x-axis ticks
    if (xticks) {
        axis(1, at = axTicks(1), labels = axTicks(1) / 1e6, 
             cex.axis = cex.axis,
             mgp = c(1.7, 0.4, 0), tcl = -0.3)
    } else if (!border) {
        axis(1, at = axTicks(1), labels = FALSE, tcl = -0.3)
    }
    
    # Add legend
    if (!is.null(legend_pos)) {
        leg <- pt.bg <- pch_leg <- title <- NULL
        
        # eQTL gene legend
        if (!is.null(eqtl_gene)) {
            leg <- levels(bg)[-1]
            pt.bg <- scheme[-1]
            pch_leg <- rep(21, length(scheme) - 1)
        } 
        # LD legend
        else if (showLD & hasLD) {
            leg <- c("0.8 - 1.0", "0.6 - 0.8", "0.4 - 0.6", 
                    "0.2 - 0.4", "0.0 - 0.2")
            title <- expression({r^2})
            pch_leg <- rep(21, 5)
            pt.bg <- rev(LD_scheme[-c(1, 7)])
        }
        
        # Add beta coefficient to legend
        if (!is.null(beta)) {
            leg <- c(leg, expression(beta > 0), expression(beta < 0))
            pch_leg <- c(pch_leg, 2, 6)
            pt.bg <- c(pt.bg, NA, NA)
        }
        
        # Add credible sets to legend
        if (has_cs) {
            # Create a data frame with CS info and p-values
            cs_info <- data.frame(
                cs_num = cs_numbers,
                genename = cs_genename,
                p_value = sapply(cs_numbers, function(cn) {
                    min(data$p[which(data$cs_num == cn)], na.rm = TRUE)
                })
            )
            
            # Sort by p-value
            cs_info <- cs_info[order(cs_info$p_value), ]
            
            # Create legend entries in sorted order
            cs_leg <- paste("CS", cs_info$cs_num, "(", cs_info$genename, ")")
            leg <- c(leg, cs_leg)
            
            # Maintain corresponding shapes order
            cs_shapes_ordered <- cs_shapes[((seq_along(cs_info$cs_num)-1) %% length(cs_shapes)) + 1]
            pch_leg <- c(pch_leg, cs_shapes_ordered)
            pt.bg <- c(pt.bg, rep(cs_color, length(cs_info$cs_num)))
        }
        
        # Draw legend
        if (!is.null(leg)) {
            legend(legend_pos, legend = leg, y.intersp = 0.96,
                   title = title,
                   pch = pch_leg,
                   pt.bg = pt.bg,
                   col = 'black',
                   bty = 'n',
                   cex = 0.8)
        }
    }
}

#' Add labels to plot
#' @noRd
add_labels <- function(lx, ly, labs, label_x, label_y, cex = 1) {
    label_x <- rep_len(label_x, length(lx))
    label_y <- rep_len(label_y, length(ly))
    dx <- diff(par("usr")[1:2]) * label_x /100
    dy <- diff(par("usr")[3:4]) * label_y /100
    dlines(lx, ly, dx, dy, xpd = NA)
    adj1 <- -sign(dx) *0.56+0.5
    adj2 <- -sign(dy) +0.5
    adj2[abs(label_x) > abs(label_y)] <- 0.5
    adj1[abs(label_x) < abs(label_y)] <- 0.5
    if (length(unique(adj1)) == 1 & length(unique(adj2)) == 1) {
        text(lx + dx, ly + dy, labs,
             adj = c(adj1[1], adj2[1]), cex = cex, xpd = NA)
    } else {
        adj <- cbind(adj1, adj2)
        for (i in seq_along(labs)) {
            text(lx[i] + dx[i], ly[i] + dy[i], labs[i],
                 adj = adj[i,], cex = cex, xpd = NA)
        }
    }
}

#' Draw connecting lines
#' @noRd
dlines <- function(x, y, dx, dy, ...) {
    mx <- cbind(x, x + dx, NA)
    my <- cbind(y, y + dy, NA)
    xs <- as.vector(t(mx))
    ys <- as.vector(t(my))
    lines(xs, ys, ...)
}
