# GPCNA: Gene Pleiotropy Co-expression Networks Analysis

## Example

The file test-data.rds contains an imaginary expression profile data set made of 20 samples and 1000 genes. Additionally, files A.txt, B.txt and C.txt contains imaginary lists of marker genes for three different types of cell. In order to analyze which genes "C" as secondary cell type lets do the following steps.

After loading the library we create a GCN from the data using createGCN.

> library(GPCNA)

> primary.net = createGCN( "test-data.rds" )

Then, we create a new expression profile removing the effects related to cell-types "A" and "B" using removePimaryEffect.

> data = removePrimaryEffect( expr.data = "test-data.rds", target.enrichment = c("A", "B"), net = primary.net )

Notice that we pass as third argument the network we have already created. Thus, removePrimaryEffect functions does not have to build it. 

Now, we create a GCN using the new data set

> secondary.net = createGCN( data )

Notice that, createGCN is flexible enough to work both with filenames and with R objects. 

Now, we use getModulesEnrichment to take a look at the enrichment of both networks

> getModulesEnrichment( primary.net )

          blue red pink      magenta       yellow green brown turquoise black
          A 1.000000e+00   1    1 5.647068e-30 1.000000e+00     1     1         1     1
          B 1.000000e+00   1    1 1.000000e+00 4.206175e-25     1     1         1     1
          C 6.402609e-20   1    1 1.000000e+00 1.000000e+00     1     1         1     1

> getModulesEnrichment( secondary.net )

          blue green red yellow brown turquoise
          A 1.000000e+00     1   1      1     1         1
          B 1.000000e+00     1   1      1     1         1
          C 2.460543e-18     1   1      1     1         1

as it can be seen the second network shows that just cell-type "C" remains as enrichment for the modules. In fact, if we count the number of genes in module blue we can see that now, there are more genes enriched by cell-type "C" than in the first network:

> table(primary.net$moduleColors)

    black      blue     brown     green   magenta      pink       red turquoise    yellow 
      109       160       104       100        71        80        92       180       104 
> table(secondary.net$moduleColors)

     blue     brown     green       red turquoise    yellow 
      186       158       130       138       261       127 

In this case we can see that the expression profile processed to remove the effect of cell-types "A" and "B" has more "C" enriched genes than the original expression profiles. Concretelly, 180-160 = 26 more genes enriched by cell-type "C". Let's see now where they come from, i.e., what cell-type had those genes in the first network.

> blue.genes = names(secondary.net$moduleColors[secondary.net$moduleColors == "blue"])

> changes = enrichmentEvolution( primary.net, secondary.net, genes = blue.genes)

blue.genes has the name of the genes included in the blue module of the second network and using enrichmentEvolution we obtain a data frame with all the information we need. Thus, to see just the genes that has "C" as secondary cell-type we do the following:

> changes[changes$primary.enrichment != changes$secondary.enrichment, ]

                 gene primary.module primary.enrichment secondary.module secondary.enrichment
          6    GEN_20          black                  -             blue                    C
          28  GEN_150          black                  -             blue                    C
          29  GEN_162          black                  -             blue                    C
          31  GEN_165          black                  -             blue                    C
          66  GEN_366      turquoise                  -             blue                    C
          75  GEN_391          brown                  -             blue                    C
          101 GEN_531      turquoise                  -             blue                    C
          124 GEN_707          black                  -             blue                    C
          128 GEN_719      turquoise                  -             blue                    C
          138 GEN_773      turquoise                  -             blue                    C
          147 GEN_856      turquoise                  -             blue                    C
          148 GEN_859      turquoise                  -             blue                    C
          161 GEN_922          black                  -             blue                    C
          162 GEN_923      turquoise                  -             blue                    C
          172  GEN_96        magenta                  A             blue                    C
          173 GEN_125        magenta                  A             blue                    C
          174 GEN_204        magenta                  A             blue                    C
          175 GEN_430        magenta                  A             blue                    C
          176 GEN_729        magenta                  A             blue                    C
          177 GEN_752        magenta                  A             blue                    C
          178 GEN_852        magenta                  A             blue                    C
          179  GEN_14         yellow                  B             blue                    C
          180 GEN_313         yellow                  B             blue                    C
          181 GEN_342         yellow                  B             blue                    C
          182 GEN_367         yellow                  B             blue                    C
          183 GEN_482         yellow                  B             blue                    C
          184 GEN_595         yellow                  B             blue                    C
          185 GEN_611         yellow                  B             blue                    C
          186 GEN_763         yellow                  B             blue                    C

We can see the 26 genes that have been found to be enriched by "C" in a second level.
