GModel
================
Rutuja Gupte

This guide is for GModel2. The main difference between GModel and
GModel2 is the GModel needs prior estimates of the proportion of
variation that is due to genetic effects. Thus, GModel2 takes
considerably longer to run. In fact, among the softwares being compared
here, this is the slowest. However, it is very easy to download and easy
to use but requires files following a highly specific format.

## Installation

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
