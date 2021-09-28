# Contributing to iqcToolbox
Thanks for taking the time to contribute!

## Did you find a bug?
- First check open and closed issues to determine if the bug is not already reported
- If no issues discuss the bug, open a new issue with the following elements.
  - Clear title and description
  - System information (OS, MATLAB version, etc.)
  - Code snippet that exposes the bug

## Do you want to write a patch for a bug?
- Fork a personal copy of the repository
- Develop a fix in your fork (be sure to follow the [coding style guidelines](https://github.com/iqcToolbox/iqcToolbox/wiki/iqcToolbox-Coding-Standards))
- Check that your fix does not break any preexisting tests (`tests.run_test_script`)
- Create a PR to the master branch, linking the issue with the reported bug

Along with reviewing code, the PR reviewer will add tests that fail without the patch but pass with the patch. Once review is complete and tests are passing, an admin will merge in your patch. Thanks!

## Do you want to write a new feature?
- Email the admin account to begin discussions on the feature. If the new feature aligns well with the toolbox, the team will request you post a new issue describing the feature.
- Fork a personal copy of the repository
- Develop the desired feature in your fork (be sure to follow [coding style guidelines](https://github.com/iqcToolbox/iqcToolbox/wiki/iqcToolbox-Coding-Standards))
- Check that your feature does not break any preexisting tests (`tests.run_test_script`)
- Create a PR to the develop branch, linking the new issue

In order to pass the PR, tests must be developed to check that the new code satisfies the desired feature (as described in the issue) AND that the new code has 100% statement coverage. Once the review is complete and tests are passing, an admin will merge in your patch.  Thanks!

## Do you want to add new tests?
- Create a new issue describing what testing is lacking and why
- Fork a personal copy of the repository
- Develop tests (be sure to follow [coding style guidelines](https://github.com/iqcToolbox/iqcToolbox/wiki/iqcToolbox-Coding-Standards))
- Create a PR to the develop branch, linking the new issue

Once the review is complete, an admin will merge in your new tests. Thanks!