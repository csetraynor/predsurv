# Create CGDS object
mycgds = cgdsr::CGDS("http://www.cbioportal.org/public-portal/")
assertthat::assert_that(all(cgdsr::test(mycgds) == "OK"))
# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = cgdsr::getCancerStudies(mycgds)[81,1]
mycaselist = cgdsr::getCaseLists(mycgds,mycancerstudy)[2,1]
# Get clinical data for the case list
clinicaldata = cgdsr::getClinicalData(mycgds,mycaselist)
# Get an easier code to read
names(clinicaldata) <- tolower(names(clinicaldata))
clinicaldata <- tibble::rownames_to_column(clinicaldata, var = "patient_id")
clinicaldata$patient_id <- gsub( "\\.","-", clinicaldata$patient_id)

rm(list = c("mycgds","mycancerstudy","mycaselist"))


devtools::use_data(clinicaldata)
