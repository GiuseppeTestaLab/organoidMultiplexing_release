
# Project

**Project Name - Long version**

Template github repository that can be used e.g. for the publication of the code when submitting a paper. 

GH Actions are set to:

1. convert the ipynb to html, create an index html page and store all the pages in a separate branch, where they can be published as static website with github pages
2. create a README index file for the repo that contains all the links to the html and ipynb files.

Feel free to fork it for "publishing" your own notebooks, and [contact me](mailto:matteo.bonfanti@fht.org) with suggestions for improvements.


An html version of the notebooks is accessible [here](https://matbonfanti.github.io/project-template/).




## Test Notebook

Links: [jupyter notebook](prova.ipynb) and [html file](https://matbonfanti.github.io/project-template/prova.html).

This is a **test notebook**. A very basic example of notebook, just to exemplify how the
GitHub Actions defined in this repository work.

If you want to add your own notebooks, just put the ipynb files in the root project folder and then
customize the "pages" section of the [config.yaml](https://github.com/matbonfanti/project-template/blob/main/resources/config.yaml)
file by including the notebooks that need to be indexed. The "notebook" attribute defines the name
of the ipynb file and the "description" attribute contains a descriptive string that will be
included in the README and html index.




## Test HTML File

Links: [html file](https://matbonfanti.github.io/project-template/prova2.html).

This is a **test HTML file**. This example shows how HTML files can be directly included in the 
final repository and linked from the repository index.



Here you can put some final remarks.
These will be appended at the end of the README / index page.


---
*Note: this README file has been generated automatically.* <br>
*Please do not modify it directly but instead work on [this config file](resources/config.yaml).*


