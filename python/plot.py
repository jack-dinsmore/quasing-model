import matplotlib.pyplot as plt
import numpy as np

betas = [0.4, 0.41000000000000003, 0.42000000000000004, 0.43, 0.44, 0.45, 0.46, 0.47000000000000003, 0.48, 0.49, 0.5]
order_params = [0.051732381184895836, 0.06700331164944556, 0.08913872831611223, 0.15970816663516466, 0.5101541498655914, 0.750468264343918, 0.814727783203125, 0.851683093655494, 0.8773756129767305, 0.8964162847047211, 0.911744640719506]
susceptibilities = [0.003807657508439915, 0.00633659670429845, 0.013305633298812373, 0.038892643426054266, 0.2958016754478537, 0.5644888724050214, 0.6589932287892988, 0.7255762007928663, 0.7734500310754263, 0.802207290485341, 0.8308364550272623]

betas = np.array(betas)
mid_betas = (betas[1:] + betas[:-1]) / 2
order_params = np.array(order_params)
susceptibilities = np.array(susceptibilities)
order_param_slope = (order_params[1:] - order_params[:-1])
susceptibility_slope = (susceptibilities[1:] - susceptibilities[:-1])

fig, (ax_param, ax_sus) = plt.subplots(figsize=(6, 8), nrows=2, sharex=True)

ax_param.scatter(betas, order_params, color="C0")
ax_param.plot(betas, order_params, color="C0")
ax_param.scatter(mid_betas, order_param_slope, edgecolor="C0", facecolor='none')
ax_param.plot(mid_betas, order_param_slope, color="C0", linestyle='dashed')
ax_param.set_ylabel("Magnetization")

ax_sus.scatter(betas, susceptibilities, color="C1")
ax_sus.plot(betas, susceptibilities, color="C1")
ax_sus.scatter(mid_betas, susceptibility_slope, edgecolor="C1", facecolor='none')
ax_sus.plot(mid_betas, susceptibility_slope, color="C1", linestyle='dashed')
ax_sus.set_xlabel("beta")
ax_sus.set_ylabel("Susceptibility")

fig.tight_layout()
fig.savefig("../figs/params.png")