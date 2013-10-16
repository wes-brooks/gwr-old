title.map = list(pindpov='Proportion of individuals in poverty', 
    logitindpov='logit of proportion of individuals in poverty',
    pag='Comparing coefficients of\nagricultural employment',
    pex='Comparing coefficients of\nmining employment',
    pman='Comparing coefficients of\nmanufacturing employment',
    pserve='Comparing coefficients of\nservices employment',
    potprof='Comparing coefficients of\nemployment in other professions',
    pwh='Comparing coefficients of\nproportion white',
    pblk='Comparing coefficients of\nproportion white',
    pind='Comparing coefficients of\nproportion indigenous',
    phisp='Comparing coefficients of\nproportion hispanic',
    metro='Comparing coefficients of\nmetro area indicator',
    pfampov='Proportion of families in poverty',
    logitfampov='logit of proportion of families in poverty',
    pfire='Comparing coefficients of employment in \nfinance, insurance, or real estate'
)


pp = list()

for (k in 1:length(predictors)) {
    v = predictors[k]
    cc = sapply(model[['GWEN']][['1970']][['model']][['models']], function(x) {x[['coef.unshrunk']][k+1]})
    cc = cbind(cc, sapply(model[['GWR']][['1970']][['model']][['models']], function(x) {x[['coef']][k+1]}))
    cc = as.data.frame(cc)
    colnames(cc) = c('GWEN', 'GWR')

    pp[[k]] = ggplot(cc) +
        aes(x=GWR, y=GWEN) +
        geom_point(size=1) + #position=position_jitter(height=0.01), 
        labs(title=title.map[[v]]) +
        theme(title=element_text(vjust=1)) +
        theme(plot.margin=unit(c(0,0,0,0),'cm'), legend.margin=unit(0,'cm')) +
        theme_bw()
}

#For the prelim slides
pdf(paste('~/git/gwr/figures/practice-talk/', yr, '-GWEN-GWR-comparison.pdf', sep=''), width=11, height=6)
brooks::multiplot(plotlist=pp, cols=3)
dev.off()

#For the prelim paper
pdf(paste('~/git/gwr/figures/poverty/', yr, '-GWEN-GWR-comparison.pdf', sep=''), width=8, height=11)
brooks::multiplot(plotlist=pp, cols=2)
dev.off()