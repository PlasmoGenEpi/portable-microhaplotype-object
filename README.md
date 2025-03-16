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
To run commands you may use good old make or the command runner [just](https://github.com/casey/just/) which is a better choice on Windows.
Use the `make` command or `duty` commands to generate project artefacts:
* `make help` or `just --list`: list all pre-defined tasks
* `make all` or `just all`: make everything
* `make deploy` or `just deploy`: deploys site
</details>

## Credits

This project was made with
[linkml-project-cookiecutter](https://github.com/linkml/linkml-project-cookiecutter).
