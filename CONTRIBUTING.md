# Contributing to OpenABF

OpenABF is maintained by [EduceLab](https://cs.uky.edu/dri) and developed in 
collaboration with our community of contributors. We welcome bug reports,
feature requests, and code contributions.

## Bug Reports and Feature Requests
If you think you have found a bug or if you would like to request a new 
feature, please check our 
[issue tracker](https://gitlab.com/educelab/OpenABF/-/issues) to make sure
that an issue has not already been opened on your topic.

## Workflow
1) Fork this repository
2) Create a new branch in the forked repository
    - If your branch addresses an issue from the Issue Tracker, please prepend 
    your branch name with the issue number (e.g. `9-fixes-a-bug`). This will 
    ensure that your branch and Merge Request are properly linked to the issue.
3) Open a Merge Request from your branch to this repository
    - If your branch is not ready to be merged, prepend `WIP:` to the Merge 
    Request title. This will keep an itinerant project manager from accidentally 
    merging an incomplete feature. When you are ready for code review, remove 
    the `WIP:` prefix. 
4) Code review
5) Merge

## Code style
This project uses `clang-format` to maintain a consistent code style 
throughout the code base. Code which has not been processed with
this tool will not be merged.

The easiest way to maintain a consistent style is to process all changes 
with `git clang-format` before they are committed:

```shell
# Stage original changes
git add .

# Run clang-format on staged changes
git clang-format 

# Stage formatting changes and commit
git add .
git commit -m "My commit message."
```

Don't worry if you forget to run this process before committing! You can 
always process an entire branch's diff by comparing against the default 
branch:

```shell
git clang-format develop
```

Unfortunately, `clang-format` does not handle all style issues. For a general 
overview of the EduceLab C++ style, please refer to our 
[C++ style guide](https://gitlab.com/educelab/style-guides/-/blob/master/C++%20Style%20Guide.md).

## License
Any changes intentionally contributed to this repository are assumed to
be licensed under the terms outlined in [LICENSE](LICENSE).

## Attribution
This project actively maintains a 
[citable record on Zenodo](https://doi.org/10.5281/zenodo.4483858). 
We are happy to list our active community of contributors as authors on 
this record. After you have contributed 10 or more commits to this project,
please open a new pull request which adds your name and ORCID to 
`.zenodo.json`.
