# Contribution guidelines for DuMu<sup>x</sup>

DuMu<sup>x</sup> and DUNE are community projects and we are happy about all forms of external contributions.
Here are many easy things that you can do, like

* read the documentation and tell us where and how it should be improved,
* try to install DuMu<sup>x</sup> on your platform and report bugs if it doesn’t work,
* fix bugs and open merge requests or send us patches.

If you decide to contribute code please read this contribution guide.

## Style guide
When contributing code to DuMu<sup>x</sup> please follow the [style guide](doc/styleguide.md). Your work will enjoy much smoother sailing if you stick to it with your changes. DuMu<sup>x</sup> is a pretty large project, and a consistent way of doing things really helps a lot when trying to find your way around the code.

## Contributing

You should get your changes to us in the following way:

* Get an account for our GitLab instance (you might need to contact us if you can't create projects at first login).
* Fork the core module that you want to contribute to, just as you would do on GitHub.
* Push your changes to your fork on some branch.
* Open a merge request using the branch you pushed your changes to as the source branch and the master of the DuMu<sup>x</sup> repository
  as the target branch. GitLab will usually ask you about opening a merge request if you browse it right after pushing to some branch.
* Follow the discussion on the merge request to see what improvements should be done to the branch before merging.
* If you have developer status you don't need to do a fork and you can create branches directly.

If you have any questions or complaints about this workflow of contributing to DuMu<sup>x</sup>, please ask on the
[DuMu<sup>x</sup> mailing list](https://listserv.uni-stuttgart.de/mailman/listinfo/dumux).

## Git
* Use Git to your advantage!
* Check out this [great tutorial](https://www.atlassian.com/git/tutorials/setting-up-a-repository) in order to learn how to use it.
* Also check out [Dune's Git best practices](https://www.dune-project.org/doc/guides/git_best_practices/).
* Everything should compile after every single commit.
* Make small commits with changes limited to a single issue / change.
* Format commit messages as follows:

    ```
    [topic] Brief description of the change

    Long description containing the status quo,
    the changes the commit introduces and why.
    ```

    where `topic` is usually a foldername, `[assembly]`, a model `[2p2c]`, or any other topic, e.g. `[cmake]`.

* Use `git rebase -i master` to update branches to the changes on the master branch.
* Feature branches should be called `feature/my-bla-feature`.
* Bugfix branches should be called `fix/issue-554`.
* Cleanup branches should be called `cleanup/remove-deprecated-bla`.
* Use lower case letters only, and hyphens to separate things.

## GitLab
* Open issues for bugs / discussions / feature requests.
* Give issues a meaningful title, focus each issue on a single topic, describe what the problem is and what you tried.
* Open merge request for changes to be merged into master or other branches. Pushing to master is disabled.
* Give merge requests a meaningful title, describe what the problem is, how this is fixed with this merge request.
* If you don't have the permissions to open branches, you might have permissions to do a fork in your own GitLab namespace and open a merge request from your fork into the dumux repository.
* Merge requests get reviewed by at least one main developer.
* If continuous integration (GitLabCI / BuildBot) is enabled for merge requests, the pipeline has to pass before merging anything to master.

## Backwards Compatibility
As a general rule, all changes added to the dumux master version should be made
such that:
*  all tests and modules using dumux will still compile, and
*  the user is warned at compile time of any interface changes.

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

``` c++
Deprecated::template coolMethod(a, b); //TODO: Replace this after release 3.n

```

There are other methods for deprecating old interfaces, please see [ (cppref/deprecated) ](https://en.cppreference.com/w/cpp/language/attributes/deprecated) for more information.

In addition, please be sure to:

*  mark all deprecated interfaces with the release after which the deprecated interface will be removed, and
*  add a detailed description of the changes made to the [changelog](CHANGELOG.md)

**In some more complicated cases,** guaranteeing backwards compatibility for all possible
cases is not feasible, or would require enormous intrusive changes.  In this case, we recommend you do the following:
1.  Organize all changes neatly in a merge request.
2.  Include a detailed description of the changes in the changelog.
3.  Within the comments section of the merge request, mark one of the core developers,
 and document the reasons why guaranteeing backwards compatibility would not be feasible, and which cases will likely be effected.
4.  The core developers will decide if the changes and efforts are sufficient during the next monthly core developers meeting (aka DumuxDay).
5.  Should your efforts be deemed sufficient, continue with the standard MR procedure.

## Patches

* Patches can be supplied via the mailing list,
* should be formatted with git format-patch.
* TODO: How to format patches
