package org.pankratzlab.supernovo.output;

import java.lang.reflect.Field;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.pankratzlab.supernovo.App;
import com.google.common.base.Optional;
import com.google.common.base.Suppliers;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.SetMultimap;

/**
 * Interface to specify and allow a class's public fields to be used as output columns instead of
 * resolving to a single column with {@link Object#toString()}
 */
public interface OutputFields {

  /** Holds static final fields that should not be printed on output */
  static class Constants {
    private Constants() {}

    static final String DELIM = "\t";
    static final String MISSING = ".";
    static final Collector<CharSequence, ?, String> JOIN_COLLECTOR = Collectors.joining(DELIM);
    static final SetMultimap<Field, String> FIELD_HEADERS = HashMultimap.create();
  }

  default String generateLine() {
    return fieldValues().collect(Constants.JOIN_COLLECTOR);
  }

  default Stream<String> fieldValues() {
    return Stream.of(this.getClass().getFields()).flatMap(this::recurseValues);
  }

  default Object getOwnField(Field field) {
    try {
      return field.get(this);
    } catch (IllegalArgumentException | IllegalAccessException e) {
      throw new IllegalStateException(e);
    }
  }

  default Stream<String> recurseValues(Field field) {
    Object value = getOwnField(field);
    if (value == null)
      return Stream.generate(Suppliers.ofInstance(Constants.MISSING))
          .limit(Constants.FIELD_HEADERS.get(field).size());
    if (value instanceof OutputFields) {
      return ((OutputFields) value).fieldValues();
    }
    if (value instanceof Optional<?>) {
      return Stream.of(((Optional<?>) value).transform(Object::toString).or(Constants.MISSING));
    }
    return Stream.of(value.toString());
  }

  default String generateHeader() {
    return fieldHeaders().collect(Constants.JOIN_COLLECTOR);
  }

  default Stream<String> fieldHeaders() {
    return Stream.of(this.getClass().getFields()).flatMap(this::recurseHeaders);
  }

  default Stream<String> recurseHeaders(Field field) {
    Object val = null;
    try {
      val = field.get(this);
    } catch (IllegalArgumentException | IllegalAccessException e) {
      App.LOG.error(e);
    }
    final Stream<String> headers;
    Class<?> fieldType = field.getType();
    if (OutputFields.class.isAssignableFrom(fieldType) && val != null) {
      headers = prefixHeaders(field, ((OutputFields) val).fieldHeaders());
    } else {
      headers = Stream.of(field.getName());
    }
    return headers.peek(header -> Constants.FIELD_HEADERS.put(field, header));
  }

  static Stream<String> prefixHeaders(Field field, Stream<? extends Object> headers) {
    final String prefix = field.getName() + "_";
    return headers.map(Object::toString).map(h -> prefix + h);
  }
}
