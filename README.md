# portable-microhaplotype-object

A schema to define the minimum amount of data needed to export a microhaplotype calling pipeline analysis with associated metadata

## Website

[https://PlasmoGenEpi.github.io/portable-microhaplotype-object](https://PlasmoGenEpi.github.io/portable-microhaplotype-object)

## Repository Structure

* [examples/](examples/) - example data
* [project/](project/) - project files (do not edit these)
* [src/](src/) - source files (edit these)
  * [portable_microhaplotype_object](src/portable_microhaplotype_object)
    * [schema](src/portable_microhaplotype_object/schema) -- LinkML schema
      (edit this)
    * [datamodel](src/portable_microhaplotype_object/datamodel) -- generated
      Python datamodel
* [tests/](tests/) - Python tests

## Developer Documentation

<details>
Use the `make` command to generate project artefacts:

* `make all`: make everything
* `make deploy`: deploys site
</details>

## Credits

This project was made with
[linkml-project-cookiecutter](https://github.com/linkml/linkml-project-cookiecutter).
