test_that("one plus one is two", {
  expect_equal(1 + 1, 2)
})

test_that("you can skip tests if needed", {
  skip("This tests hasn't been written yet")
})
