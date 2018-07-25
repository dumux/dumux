# Contribution guidelines for DuMu<sup>x</sup>

DuMu<sup>x</sup> and DUNE are community projects and we are happy about all forms of external contributions.
here are many easy things that you can do, like

* read the documentation and tell us where and how it should be improved,
* try to install DuMu<sup>x</sup> on your platform and report bugs if it doesn’t work,
* fix bugs and open merge requests or send us patches

If you decide to contribute code please read this contribution guide.

## Style guide
When contributing code to DuMu<sup>x</sup> please follow the [style guide](doc/styleguide.md). Your work will enjoy much smoother sailing if stick to it with your changes. DuMu<sup>x</sup> is a pretty large project, and a consistent way of doing things really helps a lot when trying to find your way around the code.

## Contributing

You should get your changes to us in the following way:

* Get an account for our GitLab instance (you might need to contact us if you can't create projects at first login)
* Fork the core module that you want to contribute to, just as you would do on GitHub
* Push your changes to your fork on some branch
* Open a merge request using the branch you pushed your changes to as the source branch and the master of the DuMu<sup>x</sup> repository
  as the target branch. GitLab will usually ask you about opening a merge request if you browse it right after pushing to some branch
* Follow the discussion on the merge request to see what improvements should be done to the branch before merging
* If you have developer status you don't need to do a fork and you can create branches directly

If you have any questions or complaints about this workflow of contributing to DuMu<sup>x</sup>, please ask on the
DuMu<sup>x</sup> mailing list.

## Git
* Use git to your advantage!
* Check out this [great tutorial](https://www.atlassian.com/git/tutorials/setting-up-a-repository) in order to learn how to use it
* Also check out [Dune's git best practices](https://www.dune-project.org/doc/guides/git_best_practices/)
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
* feature branches should be called `feature/my-bla-feature`
* bugfix branches should be called `fix/issue-554`
* cleanup branches should be called `cleanup/remove-deprecated-bla`
* use lower case letters only, and hyphens to separate things

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
