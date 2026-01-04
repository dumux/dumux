# Contribution guidelines

DuMu<sup>x</sup> and DUNE are community projects and we are happy about all forms of external contributions.
Here are some easy things that you can do:

* Read the documentation and tell us where and how it should be improved,
* Try to install DuMu<sup>x</sup> on your platform and report bugs if it doesn’t work,
* Fix bugs and open merge requests or send us patches.

If you decide to contribute to DuMu<sup>x</sup>, please read this contribution guide.

## Style guide

When contributing code to DuMu<sup>x</sup> please follow the [style guide](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/doc/styleguide.md). Your work will enjoy much smoother sailing if you stick to it with your changes. DuMu<sup>x</sup> is a pretty large project, and a consistent way of doing things really helps a lot when trying to find your way around the code.

## Making changes

We distinguish between those with *developer status* and *external contributors*. For these two categories, two separate workflows are defined.

### As an external contributor

As an external contributor, you need to get an account for our GitLab instance first. You might need to contact us if you cannot create projects at first login. External contributions undergo the following steps:

* Create a personal *fork* of the core module that you want to contribute to.
* Push your changes to your fork on a dedicated branch.
* Open a merge request using the branch you pushed your changes to as the source branch and the master of the DuMu<sup>x</sup> repository
  as the target branch. GitLab will usually ask you about opening a merge request if you browse it right after pushing to some branch.
* Follow the discussion on the merge request to see what improvements should be done to the branch before merging.
* Wait for your merge request to undergo a review by someone with developer status. You do not need to assign a reviewer - we check merge requests that still require a review.
* If there are any: implement suggested improvements to your contribution, and trigger the CI pipelines that you can trigger yourself.
* Wait for someone with developer status to trigger the pipelines for the lecture and the testing repository.
* If all pipelines succeed, and the reviewer approves of all your changes, your branch will be merged by a developer with merge rights.

Important note: Please make sure that at least one of the following conditions applies:

* The visibility of your personal fork is set to `public` (can be set when creating the fork, do not choose `internal` or `private`). If your existing fork's visibility is not set to `public`, you can change this under Settings->General->Visibility, project features, permissions.
* Job tokens are *allowed* for your personal fork. You can set this under Settings->CI/CD->Job token permissions->Authorized groups and projects. Here, pick the option "All groups and projects".

This is necessary for the CI pipelines to succeed. Among them are the aforementioned pipelines for the lecture and the testing repositories. These can only be run by someone with developer status.

If you have any questions or remarks about this workflow of contributing to DuMu<sup>x</sup>, please ask on the [DuMu<sup>x</sup> support channel on Matrix](https://matrix.to/#/!dKKvOPMFJwyhekAKbj:matrix.org?via=matrix.org&via=gitter.im&via=matrix.sp-codes.de) or the [DuMu<sup>x</sup> mailing list](https://listserv.uni-stuttgart.de/mailman/listinfo/dumux).

### As someone with developer status

If you have developer status, you do not need to create a fork. You can create branches directly in the corresponding core module.

## Git

* Use Git to your advantage!
* Check out this [great tutorial](https://www.atlassian.com/git/tutorials/setting-up-a-repository) in order to learn how to use it.
* Everything should compile after every single commit.
* Make small commits with changes limited to a single issue / change.
* Format commit messages as follows:

    ```
    [topic] Add brief description of the change

    Long description containing the status quo,
    the changes the commit introduces and why.
    ```

    where `topic` is usually a foldername, `[assembly]`, a model `[2p2c]`, or any other topic, e.g. `[cmake]`.
    For consistency and brevity use imperative `Add`/`Change`/`Improve`/`Remove` etc. instead of other verb forms like `Added`/`Adding`.

* Use `git rebase -i master` to update branches to the changes on the master branch.
* Feature branches should be called `feature/my-bla-feature`.
* Bugfix branches should be called `fix/issue-554`.
* Cleanup branches should be called `cleanup/remove-deprecated-bla`.
* Use lower case letters only, and hyphens to separate things.
* In case multiple authors have contributed to a single commit, list co-authors at the end of the commit message as

    ```
    [topic] Brief description

    Long description

    Co-authored-by: John Doe <john@doe.com>
    Co-authored-by: Max More <max.more@bla.org>
    ```

## GitLab

* Open issues for bugs / discussions / feature requests.
* Give issues a meaningful title, focus each issue on a single topic, describe what the problem is and what you tried.
* Open merge request for changes to be merged into master or other branches. Pushing to master is disabled.
* Give merge requests a meaningful title, describe what the problem is, how this is fixed with this merge request.
* Merge requests get reviewed by at least one main developer.
* If continuous integration (GitLabCI / BuildBot) is enabled for merge requests, the pipeline has to pass before merging anything to master.

## Backwards Compatibility

As a general rule, all changes added to the DuMu<sup>x</sup> master version should be made
such that:

* all tests and modules using DuMu<sup>x</sup> will still compile, and
* the user is warned at compile time of any interface changes.

This can be done by masking removed methods with deprecation warnings: (example 1)

```c++
[[deprecated("This method will be removed after release (3.n). Use newMethod(x,y,z) instead!")]]
int oldMethod(const int x, const int y) const
{ return x + y; }

int newMethod(const int x, const int y, const int z) const
{ return x + y > z ? x + y : z; }
```

or by adding intermediate method calls with warnings guarded by isValid
in the [`dumux/common/deprecated.hh`](dumux/common/deprecated.hh) header: (example 2)

```c++
// support old interface of the coolMethod() function on problems
template<class B>
struct HasNewCoolMethodIF
{ template<class A> auto operator()(A&& a) -> decltype( a.coolMethod(std::declval<const B>()) ) {} };

template<class A, class B
         typename std::enable_if_t<!decltype(isValid(HasNewCoolMethodIF<B>()).template check<A>())::value, int> = 0>
[[deprecated("Use new coolMethod() interface that additionally receives
              the object b! This will be removed after 3.n release")]]
auto coolMethod(const A a, const B b)
{ return a.coolMethod(); }

template<class A, class B,
         typename std::enable_if_t<decltype(isValid(HasNewCoolMethodIF<B>()).template check<A>())::value, int> = 0>
auto coolMethod(const A a, const B b)
{ return a.coolMethod(b); }
```

and replace the call with:

```c++
Deprecated::coolMethod(a, b); //TODO: Replace this after release 3.n
```

There are other methods for deprecating old interfaces, please see [(cppref/deprecated)](https://en.cppreference.com/w/cpp/language/attributes/deprecated) for more information.

In addition, please be sure to:

* mark all deprecated interfaces with the release after which the deprecated interface will be removed, and
* add a detailed description of the changes made to the [changelog](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/CHANGELOG.md)

**In some more complicated cases,** guaranteeing backwards compatibility for all possible
cases is not feasible, or would require enormous intrusive changes. In this case, we recommend you do the following:

1. Organize all changes neatly in a merge request.
2. Include a detailed description of the changes in the changelog.
3. Within the comments section of the merge request, mark one of the core developers, and document the reasons why guaranteeing backwards compatibility would not be feasible, and which cases will likely be effected.
4. The core developers will decide if the changes and efforts are sufficient during the next monthly core developers meeting (aka DumuxDay).
5. Should your efforts be deemed sufficient, continue with the standard MR procedure.

## Patches

* Patches can be supplied via the mailing list, and
* should be formatted with git format-patch.
* We recommend using GitLab to provide patches in the form of a merge request (fork workflow).
