# How to do pull requests 

You just finished your feature branch and want to merge it with the main. Please follow this procedure to make sure
that nothing has been broken:
1. Call `git merge main feature`.
2. Run all the tests. Please make sure that every single test passes.
3. Once every single test passes, you can `git checkout main` and `git merge feature main`.
4. Finally you are ready to `git push`. This will create a pull request.
5. It is helpful if you say a few words about your feature branch in the pull request. Moreover you should emphasize that you followed the procedure listed here.
6. One of us will then approve your request, assuming we don't find anything strange.
Congratulations! You are now a contributor!

