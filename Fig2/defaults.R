### defaults
N <- 1e4
S <- 1e4
S_caus <- S*0.1
h2_add <- 0.3
h2_coord <- 0 #(0.1*h2_add)/2 # shared gxe
h2_uncoord <- 0 #(0.1*h2_add)/2
h2_gxe <- 0.2
alpha <- 0
taus <- 10^seq(-10, 0, by = 1)
tuner <- 10
expected_n_iter <- 1e3 # 10 iter = 10min
eta <- 1 # heteroskedasticity
lambda_boxcox <- 1
pgs_methods <- c("pgs", "ampPGS", "PGSC")
base_cols <- c('#004488','#6699CC','#994455')
names(base_cols) <- pgs_methods
trans_cols <- c(paste0(base_cols[1],"99"),paste0(base_cols[2],"99"),
                paste0(base_cols[3],"99"))
names(trans_cols) <- pgs_methods
ltys <- c(4,3,4)
names(ltys) <- pgs_methods
