# Unit testing notes (inspired by working on getBatch)
I narrowed down what I read on the internet to what I found useful when writing unit tests for getBatch.

The idea behind unit testing is write a bunch of functions that you can run after you make changes to your code, just to check that everything is still running as it should. "Unit" in "unit tests" is often defined as "the smallest component that it makes sense to test", so it may make sense to do test for different components of a function. As such, informative test names are useful.

Typically you make tests to cover the happy paths and any edge cases you can see. A common edge case is passing in null/garbage. Usually you have 2-3 tests per method, if you find yourself writing dozens of tests per method, your method is way to complex and needs to be broken into simpler parts (most likely, obviously there are exceptions to every rule.) 

Expect errors for all cases where input checks are used to safeguard intended behavior in the unit test.

Check rather mechanically at each point of a function that the expected objects are there. After calls to helpers etc.

Make sure that you test the happy path, that is everything is working as expected including every line of a function; take care with control flow depending on arguments etc.

One of the most important feature of a unit test is that it's isolated from any other external influences. This ensures that you're only testing the code that you've written and aren't dependent on the state of other modules. In order to do this you should write test doubles to replace your production dependencies for testing purposes. These generally take the form of a stub or mock object which returns the data you specify in order to fulfill your test expectation.

Use tryCatch in places where you can't manually check arguments for the correct datatype etc.

Don't do tests about silly input (best handled by input checks) or warnings which one wants to be errors (best handled by tryCatch etc.)

https://stackoverflow.com/questions/63017378/unit-tests-and-checks-in-package-function-do-we-do-checks-in-both

Assertthat is a separate package for ascertaining the validity of parameters inside functions!

