GModel
================
Rutuja Gupte
2024-05-31

# GModel

This guide is for GModel2. The main difference between GModel and
GModel2 is the GModel needs prior estimates of the proportion of
variation that is due to genetic effects. Thus, GModel2 takes
considerably longer to run. In fact, among the softwares being compared
here, this is the slowest. However, it is very easy to download and easy
to use but requires files following a highly specific format.

The file formatting guide here is based in R using the vcfR package. I
am assuming that R is installed. If not, please visit [The Comprehensive
R Archive Network](https://cran.rstudio.com/) and find the appropriate
version of R.

## Installations

Let us begin with the installations. There will be a lot of file
processing using vcfR. I also like to use tidyverse for managing my
data. But you do not have to use tidyverse. If these libraries are not
installed. Run the following lines of code on an R console. vcfR may
take a longer time to install.

    install.packages("tidyverse")
    install.packages("vcfR")

### Troubleshooting

This will work with Windows but may not work with Linux. As a side note,
some R packages need a C and a Fortran compile to compile some packages.
I already had the GCC compiler and had to download a Fortran compiler.

    sudo apt-get install gfortran

vcfR also requires ape and vegan packages other than the standard CRAN
repositories. For Linux, if install.packages is not working. Start by
downloading the source code for all 3 packages (.tar.gz files) in the
appropriate folder with the other packages.

The source links for the packages can be found on their official CRAN
website. 1.
[ape](%22https://cran.r-project.org/web/packages/ape/index.html%22) 2.
[vegan](%22https://cran.r-project.org/web/packages/vegan/index.html%22)
3. [vcfR](%22https://cran.r-project.org/web/packages/vcfR/index.html%22)

Then install the 3 packages in this order using R code that looks like
this.

    install.packages("/home/rgupte/R/x86_64-pc-linux-gnu-library/4.4/ape_5.8.tar.gz",lib="/home/rgupte/R/x86_64-pc-linux-gnu-library/4.4", repos = NULL, type="source")
    install.packages("/home/rgupte/R/x86_64-pc-linux-gnu-library/4.4/vegan_2.6-6.1.tar.gz",lib="/home/rgupte/R/x86_64-pc-linux-gnu-library/4.4", repos = NULL, type="source")
    install.packages("/home/rgupte/R/x86_64-pc-linux-gnu-library/4.4/vcfR_1.15.0.tar.gz",lib="/home/rgupte/R/x86_64-pc-linux-gnu-library/4.4", repos = NULL, type="source")

### GModel

Now time to download GModel. Open a terminal.

    wget "https://bernardo-group.org/wp-content/uploads/2021/01/GModel.zip"
    unzip GModel.zip
    cd GModel_release/
    chmod 755 GModel2.exe

Make sure that GModel2.exe is executable and test it.

    ./Gmodel2.exe

There will be a prompt asking for the parameter file. The GModel
download folder contains sample files for both GModel and GModel2. To
check if your installation has working, pass ParmsGModel2.csv as the
parameter file. Again, warning, it can take more than an hour to run.