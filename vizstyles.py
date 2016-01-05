def periodic_color_scale(color_list):
    class PeriodicColorScale:
        """
        Color generator
        """

        COLORS = color_list

        def __init__(self, config={}):
            """
            Initialize a color scale

            Parameters
            ----------
            config : dict, default {}
                Dict for configure default colors. Its values are colors unique
                to each key. Colors included in config will never be used.

            """
            self.__config = config
            self.__buffer = self.COLORS[:]

            for color in self.__config.values():
                if color in self.__buffer:
                    self.__buffer.remove(color)

        def get_color(self, name):
            """
            Get color unique to the recieved name

            Parameters
            ----------
            name : string
                This method returns one color unique to this parameter.

            """
            if self.__config.get(name) is None:
                self.__config[name] = self.__buffer.pop(0)
                if len(self.__buffer) == 0:
                    self.__buffer = self.COLORS[:]

            return self.__config[name]

        def get_config(self):
            """Get an instance of dic as the config of colors."""
            return self.__config

    return PeriodicColorScale

elegans_color_scale = periodic_color_scale(
        ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#e31a1c", "#8dd3c7",
         "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69",
         "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f"])
matplotlib_color_scale = periodic_color_scale(
        ["#0000ff", "#008800", "#ff0000", "#ff00ff", "#ffff00", "#00ffff",
         "#000000"])
default_color_scale = elegans_color_scale
