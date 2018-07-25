# Contribution guidelines for DuMu<sup>x</sup>

## Style guide
When contributing code to DuMu<sup>x</sup> please follow the [styleguide](doc/styleguide.md).

## Git
* Everything should compile after every single commit
* Make small commits with changes limited to a single issue / change
* Format commit messages as follows
```
[topic] Brief description of the change

Long description containing the status quo,
the changes the commit introduces and why.
```
where `topic` is usually a foldername, `[assembly]`, a model `[2p2c]`, or any other topic, e.g. `[cmake]`.

* Use `git rebase -i master` to update branches to the changes on the master branch

## GitLab
* open issues for bugs / discussions / feature requests
* give issues a meaningful title, focus each issue on a single topic, describe what the problem is and what you tried
* open merge request for changes to be merged into master or other branches. Pushing to master is disabled.
* give merge requests a meaningful title, describe what the problem is, how this is fixed with this merge request
* if you don't have the permissions to open branches, you might have permissions to do a fork in your own GitLab namespace and open a merge request from your fork into the dumux repository.
* merge requests get reviewed by at least one main developer
* if continuous integration (GitLabCI / BuildBot) is enabled for merge requests, the pipeline has to pass before merging anything to master

## Patches

* Patches can be supplied via the mailing list
* Should be formatted with git format-patch
* TODO: How to format patches
