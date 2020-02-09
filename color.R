
## Henry Linder 2020 mhlinder@gmail.com

library(magrittr)
library(dplyr)

df_cols <-
    c(
        ## Turquoise green
        "#1abc9c", "Turquoise",     "GreenA1",
        "#16a085", "Green sea",     "GreenA2",
        ## Green
        "#2ecc71", "Emerald",       "GreenB1",
        "#27ae60", "Nephritis",     "GreenB2",
        ## Blue
        "#3498db", "Peter river",   "Blue1",
        "#2980b9", "Belize hole",   "Blue2",
        ## Purple
        "#9b59b6", "Amethyst",      "Purple1",
        "#8e44ad", "Wisteria",      "Purple2",
        ## "Black"
        "#34495e", "Wet asphalt",   "Black1",
        "#2c3e50", "Midnight blue", "Black2",
        ## Yellow
        "#f1c40f", "Sun flower",    "Yellow1",
        "#f39c12", "Orange",        "Yellow2",
        ## Orange
        "#e67e22", "Carrot",        "Orange1",
        "#d35400", "Pumpkin",       "Orange2",
        ## Red
        "#e74c3c", "Alizarin",      "Red1",
        "#c0392b", "Pomegranate",   "Red2",
        ## White, grey 1, grey 2, grey 3
        "#ecf0f1", "Clouds",        "White",
        "#bdc3c7", "Silver",        "Grey1",
        "#95a5a6", "Concrete",      "Grey2",
        "#7f8c8d", "Asbestos",      "Grey3"
    ) %>%
    matrix(ncol = 3, byrow = TRUE) %>%
    as.data.frame(stringsAsFactors = FALSE)
names(df_cols) <- c("Hex", "Name", "ID")
##
cols_all <- df_cols$Hex
names(cols_all) <- df_cols$ID
##
reordered <- c("Red1", "Orange1", "Yellow1", "GreenB1", "GreenA1", "Blue1", "Purple1", "Black1",
               "Red2", "Orange2", "Yellow2", "GreenB2", "GreenA2", "Blue2", "Purple2", "Black2",
               "Grey3", "Grey2", "Grey1", "White")
cols_all <- cols_all[reordered]

pad_colors <- function(n, cols = cols_all) {
    k <- floor(n / length(cols))
    l <- n %% length(cols)
    cols <- c(rep(cols, k), cols[1:l])
    cols
}
