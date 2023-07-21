def image_restoration(n,m):

        n = 5000
        m = 1000

        minus1 = ConstantScalar(-1.0)
        lambda_ = Scalar("lambda")
        lambda_.set_property(Property.POSITIVE)
        sigma_ = Scalar("sigma_sq")
        sigma_.set_property(Property.POSITIVE)

        H = Matrix("H", (m, n), properties = [Property.FULL_RANK])
        I = IdentityMatrix(n, n)

        v_k = Vector("v_k", (n, 1))
        u_k = Vector("u_k", (n, 1))
        y = Vector("y", (m, 1))
        x = Vector("x", (n, 1))


        equations = Equations(
                            Equal(
                                x,
                                Times(
                                    # (H^t * H + lambda * sigma^2 * I_n)^-1
                                    Inverse( Plus(
                                        Times(
                                            Transpose(H),
                                            H
                                        ),
                                        Times(
                                            lambda_,
                                            sigma_,
                                            I
                                        )
                                    )),
                                    # (H^T * y + lambda * sigma^2 * (v - u))
                                    Plus(
                                        Times(
                                            Transpose(H),
                                            y
                                        ),
                                        Times(
                                            lambda_,
                                            sigma_,
                                            Plus(
                                                v_k,
                                                Times(
                                                    minus1,
                                                    u_k
                                                )
                                            )
                                        )
                                    )
                                )
                            )
                        )

        return equations