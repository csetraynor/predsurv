# Package predsurv

This package includes relevant code for building relevant prognostic models in medicine.
The most usual modeling techniques are available as well as measures to asses model performance.
The main functions are wrapped in module blocks mostly from relevant packages such as prodlim, pec, and resample.
This package is perfect for visualising the relevant measures.

PS:It is planned to develop also measures for model performance with time-dependent covariates. 

## Getting Started

Before installing the package make sure that you have all the dependencies shown at next section.
I would not recommend install the package without installing first all the dependencies because it can take long to install specially 
devtools::install_github("csetraynor/predsurv")

### Prerequisites

This package is build in R to install R follow R-Cran guidlines.


Many functions wraps many important packages for survival analysis.
For learning more about survival analysis in R: R-Cran Task survival

```
		snowball,
		purrr,
		caret,
		dplyr,
		ggplot2,
		survival,
		glmnet,
		pec
		
```

### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc

