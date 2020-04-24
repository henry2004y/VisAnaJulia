# Text Data Processing

Now it is an age of big data.
The fundamental of big data is how to process data efficiently.
This section is about text data processing, including data organization, regular expression, and web scraping.

Examples can be found in the [timing](https://github.com/henry2004y/VisAnaJulia/tree/master/timing) folder.

## Data Format

There are two basic packages for reading/writing ascii data:
* [DelimitedFiles](https://docs.julialang.org/en/v1/stdlib/DelimitedFiles/)
* [CSV](https://juliadata.github.io/CSV.jl/stable/)

The first one is a built-in package for handling relatively simple data, while the second one is for handling complicated data with seamless support for [DataFrame](https://juliadata.github.io/DataFrames.jl/stable/).

## Regular Expression

Julia has built-in support for regular expressions.
The best reference and tool I have ever used for regular expression in general is [RegExr](https://regexr.com/).

## Web Scraping

Packages worth mentioning:
* [HTTP](https://github.com/JuliaWeb/HTTP.jl)
* [Gumbo](https://github.com/JuliaWeb/Gumbo.jl)
* [Cascadia](https://github.com/Algocircle/Cascadia.jl)
